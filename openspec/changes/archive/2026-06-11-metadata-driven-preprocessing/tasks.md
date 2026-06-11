### Task 1: Create node/scripts/pipeline.R

**Goal:** New module with `detect_pipeline()`, `is_quantile_normalized()`, `apply_pipeline_transform()`

**Files:**
- Create: `node/scripts/pipeline.R`
- Create: `tests/testthat/test-pipeline.R`

---

### Task 2: Update node/scripts/normalize.R

**Goal:** Replace `normalize_expr_matrix()` with pipeline-driven delegation

**Files:**
- Modify: `node/scripts/normalize.R`
- Modify: `tests/testthat/test-normalize.R`

---

### Task 3: Update node/scripts/fetch.R

**Goal:** Call `detect_pipeline()` after pData extraction, pass to `apply_pipeline_transform()`

**Files:**
- Modify: `node/scripts/fetch.R`

---

### Task 4: Update node/scripts/qc.R

**Goal:** Report QN status from `is_quantile_normalized()`

**Files:**
- Modify: `node/scripts/qc.R`

---

### Task 5: Full test suite verification

**Goal:** All tests pass, end-to-end smoke test on GSE318047 + GSE100155
