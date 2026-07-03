# Metabolic Signature Discovery Platform — Preview

แพลตฟอร์มค้นหา "ลายพิมพ์ชีวภาพ" (metabolic signature) จากข้อมูล NMR metabolomics
เว็บหน้าเดียว (`index.html`) + Backend Python ที่รันไพป์ไลน์วิเคราะห์ 8 ชั้น

---

## ภาพรวม

```
┌────────────────┐      fetch /api/metabolites     ┌──────────────────────┐
│  index.html    │ ─────────────────────────────►  │  FastAPI (port 8080) │
│  (frontend UI) │ ◄─────  ผลวิเคราะห์จริง (JSON)  ──── │  api/main.py         │
└────────────────┘                                  └──────────┬───────────┘
        │ ถ้าไม่มี backend → ใช้ mock ในไฟล์เอง                    │
        │ (ยังเปิดดูได้แบบ standalone)                    run_pipeline()
                                                               ▼
                                              Layer 0 → 1 → 2 → … → 8
```

- **Frontend**: HTML/CSS/JS ล้วน (Tailwind CDN) — 9 แท็บ: ภาพรวม, จับคู่รหัสสาร, สถิติ, เครือข่าย, คลัสเตอร์, เส้นทางชีวเคมี, ปัญญาประดิษฐ์/SHAP, แผนที่ลายพิมพ์, ส่งออก
- **Backend**: Python (numpy/scipy/pandas/FastAPI) — คำนวณจริงแล้วส่งค่าให้ frontend

---

## Pipeline 8 ชั้น

| Layer | ทำอะไร | สถานะ |
|-------|--------|-------|
| 0 Ingestion | อ่านข้อมูล, melt, หา orphan | ✅ จริง |
| 1 Resolution | map รหัสสาร → `chebi_id` + **Alanine gate check** | ✅ จริง |
| 2 Normalization | PQN → log2(x+1) → z-score + **QC checkpoint** | ✅ จริง |
| 3 Differential | Shapiro→t-test/Mann-Whitney, Cohen's d, BH-FDR | ✅ จริง |
| 4 Network | correlation (Pearson+Spearman), degree, hub | ◆ deterministic |
| 5 Modules | แยกคลัสเตอร์ (spectral) → M1/M2 | ◆ deterministic |
| 6 Pathway | KEGG enrichment (hypergeometric) | ◆ deterministic |
| 7 ML / SHAP | ความสำคัญ feature (AUC proxy) | ◆ deterministic |
| 8 Evidence | รวมคะแนน → convergence score + tier | ◆ deterministic |

> ◆ = ทำงานได้จริง output ตรง shape ของจริง แต่ยังไม่ใช้ XGBoost/Louvain ตัวเต็ม (เติมทีหลังได้)

---

## ผลลัพธ์ตัวอย่าง (ข้อมูลสังเคราะห์ปัจจุบัน)

- **Gate + QC**: ผ่าน (Alanine ครบ 5 แหล่ง, z-median ≈ −0.14)
- **Modules**: M1 Energy (5 สาร) / M2 Nitrogen Shunt (4 สาร) — ARI 1.0
- **Top signatures**: L-Glutamine 0.93, L-Glutamic acid 0.90, L-Lactic acid 0.80 (High Tier)
- **Pathway**: D-Glutamine/glutamate metabolism enriched (p=0.083)

---

## วิธีรัน

```bash
cd backend
pip install -r requirements.txt
python -m uvicorn api.main:app --port 8080
# เปิด http://127.0.0.1:8080/
```

รันไพป์ไลน์เดี่ยว ๆ ดู JSON: `python pipeline.py`

**API**: `/api/health` · `/api/pipeline` · `/api/metabolites` · `/api/modules` · `/api/pathways` · `POST /api/pipeline/rerun`

---

## "ข้อมูลจริงทีหลัง" — ทำแค่ 3 ขั้น

ทุก layer คุยกันผ่าน `chebi_id` (schema กลาง) → สลับแหล่งข้อมูลได้โดยไม่แตะ layer 1-8

1. เขียน `ingestion/nmr_files.py :: NmrFileDataSource.load()` ให้อ่านไฟล์จริง
2. `pipeline.py` → `run_pipeline(source=NmrFileDataSource("data/"))`
3. เติม synonym ใน `reference/metabolite_reference.json` ถ้าชื่อไม่ตรง

---

## อัปเกรดเป็น full pipeline (ทีหลัง)

uncomment deps ใน `requirements.txt` แล้วสลับทีละ layer:
`statsmodels` (mixed-effects) · `python-louvain` (Louvain จริง) · `gseapy` (KEGG) · `xgboost`+`shap`+`optuna` (nested-CV + TreeSHAP)

รายละเอียดเต็ม: [backend/README.md](backend/README.md)
