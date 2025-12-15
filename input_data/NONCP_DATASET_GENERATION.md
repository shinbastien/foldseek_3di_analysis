# Non-CP Homolog Pairs Dataset Generation

## Overview
`noncp_homolog_pairs.tsv`는 Circular Permutation (CP) 검출 벤치마크를 위한 고품질 비-CP 호몰로그 쌍 데이터셋입니다. 이 문서는 데이터셋 재생성 절차를 설명합니다.

## Prerequisites
- **필수 파일**:
  - `scope40_domains.tsv`: SCOPe 40% 도메인 메타데이터
  - `input_data/datasets/cp_positive_pairs.tsv`: CP 양성 쌍 목록 (제외용)
  - `scope_pdb/`: PDB 도메인 파일 디렉토리
- **필수 소프트웨어**:
  - Foldseek (구조 검색)
  - TM-align (구조 정렬 및 CP 검출)
  - Python 3.8+ (pandas, numpy)

## Generation Pipeline

### Step 1: Query Selection
```bash
# 기본 설정으로 2000개 쿼리 도메인 선택
# - Superfamily당 균등 샘플링
# - CP 양성 쌍 제외
python3 build_noncp_queries.py
```

**출력**:
- `noncp_queries/noncp_query_ids.txt`
- `noncp_queries/noncp_db_ids.txt`
- `noncp_queries/noncp_query_paths.txt`
- `noncp_queries/noncp_db_paths.txt`
- `noncp_queries/sf_stats.tsv`

### Step 2: Foldseek Search
```bash
# Foldseek 구조 검색 실행 (쿼리당 K=2 히트)
cd noncp_work
foldseek easy-search \
  ../noncp_queries/noncp_query_paths.txt \
  ../noncp_queries/noncp_db_paths.txt \
  foldseek_out.0 tmp \
  --format-output query,target,fident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits
```

### Step 3: TM-align Filtering
```bash
# TM-align로 CP 시그널 검증 및 필터링
python3 scripts/run_tmalign_and_filter_noncp.py \
  --candidate_pairs noncp_work/candidate_noncp_pairs.tsv \
  --pdb_dir scope_pdb \
  --output_with_tm noncp_work/noncp_with_tm.tsv \
  --output_final input_data/noncp_homolog_pairs.tsv \
  --tm_threshold 0.5 \
  --delta_cp_threshold 0.02 \
  --n_final 2000 \
  --max_pairs_per_query 2
```

### Full Pipeline (Automated)
```bash
# 전체 파이프라인 실행 (권장)
bash scripts/run_noncp_pipeline.sh \
  --n_query 2000 \
  --n_final 2000 \
  --tm_threshold 0.5 \
  --delta_cp_threshold 0.02 \
  --k_hits 2 \
  --work_dir noncp_work
```

## Filtering Criteria

### 1. Structural Similarity
- **TM_nocp ≥ 0.5**: 명확한 구조적 호몰로그
- TM-align 정상 모드로 측정

### 2. Non-CP Verification
- **|ΔCP| ≤ 0.02**: CP 시그널 부재 확인
- ΔCP = TM_cp - TM_nocp (TM-align `-cp` 옵션 사용)
- CP가 있다면 TM_cp > TM_nocp이므로 ΔCP가 큼

### 3. Diversity
- 쿼리당 최대 2개 타겟 (중복 방지)
- Superfamily 균등 샘플링
- CP 양성 쌍 명시적 제외

## Output Format
`noncp_homolog_pairs.tsv` 컬럼:
- `query_id`, `target_id`: SCOPe 도메인 ID
- `TM_nocp`, `TM_cp`: TM-align 점수 (정상/CP 모드)
- `ΔCP`: TM_cp - TM_nocp
- `foldseek_bitscore`, `foldseek_evalue`: Foldseek 메트릭
- `query_sccs`, `target_sccs`: SCOP 분류
- `query_species`, `target_species`: PDB ID
- `query_length`, `target_length`: 서열 길이

## Current Dataset Statistics
- **Total pairs**: 1,593
- **Unique queries**: 847
- **Unique targets**: 1,521
- **Mean TM_nocp**: 0.76
- **Mean |ΔCP|**: 0.003
- **Superfamilies**: 309

## Cache and Intermediate Files
생성 과정에서 생성되는 캐시 디렉토리:
- `noncp_work/`: Foldseek 결과 및 중간 TSV (386MB)
  - `foldseek_out.*`: 병렬 Foldseek 검색 결과
  - `candidate_noncp_pairs.tsv`: 필터링 전 후보
  - `noncp_with_tm.tsv`: TM-align 결과 포함 전체 목록
- `noncp_work/tmp/`: Foldseek 임시 파일

**재생성 후 캐시 삭제 가능**: 최종 `noncp_homolog_pairs.tsv`만 유지하면 됩니다.

## Verification
재생성한 데이터셋 검증:
```bash
# 행 수 확인 (헤더 제외)
wc -l input_data/noncp_homolog_pairs.tsv
# Expected: 1594 (1593 pairs + 1 header)

# TM_nocp 및 ΔCP 통계
python3 -c "
import pandas as pd
df = pd.read_csv('input_data/noncp_homolog_pairs.tsv', sep='\t')
print(f'Pairs: {len(df)}')
print(f'Mean TM_nocp: {df[\"TM_nocp\"].mean():.4f}')
print(f'Mean |ΔCP|: {df[\"ΔCP\"].abs().mean():.6f}')
print(f'Max |ΔCP|: {df[\"ΔCP\"].abs().max():.6f}')
"
```

## Notes
- 데이터셋 생성은 TM-align 실행으로 인해 시간 소모적 (~수 시간)
- 병렬화 옵션은 `run_noncp_pipeline.sh`에서 조정 가능
- CP 양성 쌍은 `input_data/datasets/cp_positive_pairs.tsv`에서 관리

---

**Last Updated**: 2025-12-15  
**Pipeline Version**: 1.0
