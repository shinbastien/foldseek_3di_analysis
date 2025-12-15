# Input Data

## 구성

```
input_data/
├── NONCP_DATASET_GENERATION.md   # Non-CP 데이터셋 생성 문서
├── README.md                      # 이 파일
└── datasets/                      # 입력 데이터셋
    └── noncp_homolog_pairs.tsv   # Non-CP 도메인 쌍 (1593개)
```

## noncp_homolog_pairs.tsv

**위치**: `datasets/noncp_homolog_pairs.tsv`

**형식**: TSV (Tab-Separated Values)

**컬럼**:
- `query_id`, `target_id`: SCOPe 도메인 ID
- `TM_nocp`, `TM_cp`: TM-align 점수 (정상/CP 모드)
- `ΔCP`: TM_cp - TM_nocp (CP 시그널)
- `foldseek_bitscore`, `foldseek_evalue`: Foldseek 메트릭
- `query_sccs`, `target_sccs`: SCOP 분류
- `query_species`, `target_species`: PDB ID
- `query_length`, `target_length`: 서열 길이

**통계**:
- 행 수: 1593 쌍
- Unique queries: 847
- Unique targets: 1,521
- Mean TM_nocp: 0.76
- Mean |ΔCP|: 0.003 (CP 신호 없음)

**설명**: 고품질 비-CP 호몰로그 쌍 세트
- TM_nocp ≥ 0.5 필터링 (구조적 유사도)
- |ΔCP| ≤ 0.02 필터링 (CP 신호 부재)

## 재생성

`NONCP_DATASET_GENERATION.md` 참조. 완전히 재생성 가능 (deterministic pipeline, `random_state=42` 고정).

## 사용처

- `experiments/03_x1_x2_cp_detection/`: X1/X2 정렬 비교 분석
- `experiments/02_8f_10f_comparison/`: 모델 비교 벤치마크
