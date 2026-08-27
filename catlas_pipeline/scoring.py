import gc

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F


def score_with_nt(
    seq_df: pd.DataFrame,
    out_scores_path: str,
    out_embeddings_path: str,
    model_id: str = "InstaDeepAI/nucleotide-transformer-2.5b-multi-species",
    batch_size: int = 8,
    max_len: int = 512,
) -> pd.DataFrame:
    from transformers import AutoModel, AutoTokenizer

    device = "cuda" if torch.cuda.is_available() else "cpu"
    tokenizer = AutoTokenizer.from_pretrained(model_id, trust_remote_code=True)
    model = (
        AutoModel.from_pretrained(model_id, trust_remote_code=True).to(device).eval()
    )
    print(f"NT-2.5B loaded on {device}")

    def embed(seq_list, model):
        inputs = tokenizer(
            seq_list,
            return_tensors="pt",
            padding=True,
            truncation=True,
            max_length=max_len,
        )
        inputs = {k: v.to(device) for k, v in inputs.items()}
        with torch.no_grad():
            outputs = model(**inputs)
        mask = inputs["attention_mask"].unsqueeze(-1).float()
        return (outputs.last_hidden_state * mask).sum(1) / mask.sum(1)

    if not torch.cuda.is_available():
        seq_df = seq_df.head(2)  # CPU smoke-test only

    n = len(seq_df)
    embeddings_all, records = [], []
    for i in range(0, n, batch_size):
        batch = seq_df.iloc[i : i + batch_size]
        emb_ref = embed(batch["sequence_ref"].tolist(), model)
        emb_alt = embed(batch["sequence_alt"].tolist(), model)
        sim = F.cosine_similarity(emb_ref, emb_alt, dim=1)
        score = -(1 - sim).cpu().numpy()
        embeddings_all.append(emb_ref.cpu().numpy())
        for j, (_, row) in enumerate(batch.iterrows()):
            records.append(
                {
                    "variant": row["variant"],
                    "chr": row["chr"],
                    "pos": row["pos"],
                    "pval": row["pval"],
                    "nt_score": float(score[j]),
                }
            )
        if i % (batch_size * 25) == 0:
            print(f"  {i:,}/{n:,} scored...", flush=True)

    del model
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    embeddings = np.vstack(embeddings_all)
    np.save(out_embeddings_path, embeddings)

    scores_df = pd.DataFrame(records)

    scores_df = scores_df.drop_duplicates(subset="variant", keep="last")
    scores_df.to_csv(out_scores_path, sep="\t", index=False)
    print(f"Saved -> {out_scores_path}")
    return scores_df


def score_with_evo2(
    seq_df: pd.DataFrame, out_scores_path: str, model_name: str = "evo2_7b"
) -> pd.DataFrame:
    from evo2 import Evo2

    model = Evo2(model_name)

    if not torch.cuda.is_available():
        seq_df = seq_df.head(2)  # CPU smoke-test only

    records = []
    n = len(seq_df)
    for i, row in seq_df.iterrows():
        score_ref = model.score_sequences([row["sequence_ref"]])[0]
        score_alt = model.score_sequences([row["sequence_alt"]])[0]
        records.append(
            {
                "variant": row["variant"],
                "chr": row["chr"],
                "pos": row["pos"],
                "pval": row["pval"],
                "evo2_delta": float(score_alt - score_ref),
            }
        )
        if i % 25 == 0:
            print(f"  {i:,}/{n:,} scored...", flush=True)

    del model
    gc.collect()
    if torch.cuda.is_available():
        torch.cuda.empty_cache()

    scores_df = pd.DataFrame(records)
    scores_df = scores_df.drop_duplicates(subset="variant", keep="last")
    scores_df.to_csv(out_scores_path, sep="\t", index=False)
    print(f"Saved -> {out_scores_path}")
    return scores_df
