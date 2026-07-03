# Backend — Metabolic Signature Discovery Platform

Backend หลังบ้านของ `index.html` สร้างตามสถาปัตยกรรมใน Devflow doc
รันไพป์ไลน์ 8 layer แล้ว expose ผลให้ frontend ผ่าน REST API

## สรุปสั้น

- **Layer 0-3 รันจริง**: Ingestion → Resolution (+ Alanine gate) → Normalization
  (PQN + log2 + z-score + QC) → Differential (Shapiro→t/Mann-Whitney, Cohen's d, BH-FDR)
- **Layer 4-8 deterministic**: Network (Pearson+Spearman) → Modules (spectral bisection)
  → Pathway (hypergeometric) → ML importance (AUC proxy) → Evidence integration
  (min-max + equal-weight convergence score + tier cap)
- **ข้อมูลสังเคราะห์ก่อน, map ไฟล์จริงทีหลัง** ผ่าน `ingestion/` adapter — ไม่ต้องแตะ layer อื่น

## รัน

```bash
cd backend
pip install -r requirements.txt          # numpy, scipy, pandas, fastapi, uvicorn

# รันไพป์ไลน์ครั้งเดียว ดูผลเป็น JSON
python pipeline.py

# หรือรัน API server (เสิร์ฟ index.html + REST ที่ port 8080)
python -m uvicorn api.main:app --host 127.0.0.1 --port 8080 --reload
```

เปิด http://127.0.0.1:8080/ — frontend จะ `fetch('/api/metabolites')` เอาผลจริงมาแสดง
ถ้าเปิด `index.html` ตรงๆ โดยไม่มี backend จะ fallback ไปใช้ mock ในไฟล์เอง (ยังทำงานได้)

## API

| Endpoint | ผล |
|---|---|
| `GET /` | เสิร์ฟ index.html |
| `GET /api/health` | สถานะ + gate/QC ผ่านไหม |
| `GET /api/pipeline` | ผลเต็ม (metabolites, modules, pathways, log) |
| `GET /api/metabolites` | เฉพาะ metData ที่ frontend ใช้ render |
| `GET /api/modules` | การ์ด module M1/M2 |
| `GET /api/pathways` | ตาราง pathway enrichment |
| `POST /api/pipeline/rerun` | ล้าง cache แล้วคำนวณใหม่ |

## โครงสร้าง

```
backend/
  schema.py                 # data contract กลาง (ล็อกก่อน = หัวใจของ "map ทีหลัง")
  pipeline.py               # orchestrator เรียง Layer 0-8
  reference/
    metabolite_reference.json   # local CHEBI/KEGG cache (Layer 1 resolve ที่นี่)
  ingestion/
    base.py                 # DataSource interface
    synthetic.py            # ← ใช้อยู่ตอนนี้ (สร้างข้อมูลสังเคราะห์ reproducible)
    nmr_files.py            # ← STUB เติมทีหลังเมื่อมีไฟล์จริง
  layers/
    l1_resolution.py … l8_evidence.py
  api/main.py               # FastAPI
```

## วิธี "map ข้อมูลจริง" ทีหลัง (สิ่งเดียวที่ต้องทำ)

1. เขียน `ingestion/nmr_files.py::NmrFileDataSource.load()` ให้อ่านไฟล์จริง
   คืนค่าเป็น `RawCohort` (long form: `feature_label, sample_id, value, group`)
2. ใน `pipeline.py` เปลี่ยน:
   ```python
   from ingestion.nmr_files import NmrFileDataSource
   out = run_pipeline(source=NmrFileDataSource("data/"))
   ```
3. เพิ่ม synonym ใน `reference/metabolite_reference.json` ถ้า label ไฟล์จริงไม่ตรง

Layer 1-8 ไม่ต้องแก้เลย เพราะทุก layer คุยกันผ่าน `chebi_id` ตาม schema

## อัปเกรดเป็น "full pipeline" (ทีหลัง)

Layer 4-8 ตอนนี้เป็น deterministic stand-in ที่ output shape ตรงกับของจริง
สลับได้ทีละตัวโดยไม่กระทบ frontend — uncomment deps ใน `requirements.txt`:
- Layer 3 Obesity: `statsmodels` MixedLM (mixed-effects)
- Layer 5: `python-louvain` (Louvain จริง + resolution sweep)
- Layer 6: `gseapy` (KEGG enrichment)
- Layer 7: `xgboost` + `shap` + `optuna` (nested CV + TreeSHAP)
