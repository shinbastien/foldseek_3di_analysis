# Foldseek 3Di Analysis Project - 프로젝트 구조

## 📂 **전체 구조**

```
foldseek_new_3di/
├── core/                            # 핵심 인프라 & 데이터
│   ├── encoders_and_tools/          # 모델 & 도구
│   │   ├── training_3di_gpu_8f/     # 8f 모델 (encoder.pt, decoder.pt)
│   │   ├── training_3di_gpu_9f/     # 9f 모델
│   │   ├── training_3di_gpu_10f/    # 10f 모델
│   │   └── ssw/                     # SSW 바이너리 & 매트릭스
│   ├── scope_databases/             # 도메인 데이터베이스
│   │   ├── scope_pdb/               # PDB 도메인 구조 (13,712 files)
│   │   ├── scope40/                 # SCOPe 40% 비중복 데이터
│   │   └── cirpin/                  # CP 데이터베이스
│   └── foldseek_tools/              # Foldseek 유틸리티
│
├── datasets/                        # 입력 데이터셋 (중앙 저장소)
│   ├── cp_positive_pairs.tsv        # CP 양성 쌍 (451개)
│   ├── noncp_homolog_pairs.tsv      # 비-CP 호몰로그 (1,593개)
│   ├── scope40_domains.tsv          # 도메인 메타데이터
│   └── datasets/                    # 추가 파일들
│
├── experiments/                     # 실험 프로젝트들
│   ├── 01_encoder_training/         # 모델 학습 (완료)
│   │   ├── README.md
│   │   ├── results/
│   │   │   ├── training_logs_8f.txt
│   │   │   ├── training_logs_9f.txt
│   │   │   └── training_logs_10f.txt
│   │   └── notebooks/               # 학습 과정 분석
│   │
│   ├── 02_8f_10f_ssw_comparison/    # 8f vs 10f SSW 비교 (완료)
│   │   ├── README.md
│   │   ├── scripts/
│   │   │   ├── compare_8f_10f_ssw_v2.py
│   │   │   └── run_ssw_comparison.sh
│   │   ├── results/
│   │   │   ├── ssw_comparison_8f_vs_10f.csv (267MB)
│   │   │   ├── coverage95_8f_vs_10f.png
│   │   │   └── ssw_comparison_8f_vs_10f.png
│   │   └── logs/
│   │
│   ├── 03_x1_x2_cp_detection/       # X1/X2 정렬로 CP 감지 (완료)
│   │   ├── README.md
│   │   ├── scripts/
│   │   │   ├── generate_x1_x2_with_logging.py
│   │   │   ├── noncp_x1x2_compare.py
│   │   │   └── noncp_8f_10f_x1x2_compare.py
│   │   ├── check_cp/                # CP 양성 쌍 분석
│   │   │   ├── sequences/           # 3Di 시퀀스 (8f, 9f, 10f)
│   │   │   ├── results/
│   │   │   │   └── comparison_results.csv
│   │   │   └── plots/
│   │   ├── check_noncp/             # 비-CP 호몰로그 분석
│   │   │   ├── sequences/           # 3Di 시퀀스 (8f, 9f, 10f)
│   │   │   ├── results/
│   │   │   │   └── noncp_8f_10f_comparison.csv
│   │   │   └── plots/
│   │   └── results/
│   │       └── noncp_8f_10f_comparison.csv
│   │
│   └── 04_pdz_sanity_check/         # PDZ CP 검증 (완료)
│       ├── README.md                # 상세 실험 설명
│       ├── RESULTS_SUMMARY.md       # 결과 요약 & 해석
│       ├── 1_raw_3di_sequences/     # 원본 3Di 시퀀스
│       ├── 2_untrimmed_ssw_comparison/
│       ├── 3_structural_alignment/  # TM-align 결과
│       ├── 4_trimmed_ssw_comparison/
│       ├── pdb_structures/
│       ├── log.txt                  # 파이프라인 로그
│       └── *.mat                    # SSW 매트릭스
│
├── scripts/                         # 공용 스크립트
│   ├── pairwise_3di_pipeline.py     # 메인 파이프라인
│   ├── pdb_to_3di.py                # PDB → 3Di 변환
│   ├── tmalign_3di_match_pipeline.py
│   ├── ssw_identify.py
│   └── ...
│
├── utils/                           # 헬퍼 함수
│   ├── __init__.py
│   └── ...
│
├── analysis/                        # 범용 분석 결과
│   └── cp_set_analysis.tsv
│
├── docs/                            # 문서
│   ├── PROJECT_INDEX.md             # 이 파일
│   ├── README.md
│   └── SETUP.md
│
├── input_data/                      # 레거시 입력 데이터 (DEPRECATED)
│   └── datasets/
│       └── noncp_homolog_pairs.tsv
│
└── tmp/                             # 임시 작업 디렉토리
    ├── *_vs_*_w*/                   # 개별 실험 임시 파일
    └── pairwise_3di_*/
```

---

## 🗂️ **주요 디렉토리 설명**

### **core/** - 핵심 인프라

#### `encoders_and_tools/`
- **용도**: 학습된 3Di 인코더 및 필수 도구 저장
- **내용**:
  - `training_3di_gpu_8f/`: encoder.pt, decoder.pt, states.txt
  - `training_3di_gpu_9f/`: 동일 구조
  - `training_3di_gpu_10f/`: 동일 구조
  - `ssw/`: SSW 바이너리, 매트릭스 (s_8f.mat, s_10f.mat, sub_score.mat)
- **중요**: 모든 실험에서 재사용됨, **삭제 금지**

#### `scope_databases/`
- **용도**: 도메인 구조 및 메타데이터
- **내용**:
  - `scope_pdb/`: 13,712개 PDB 도메인 파일 (1.5G)
  - `scope40/`: SCOPe 40% 비중복 데이터 (638M)
  - `cirpin/`: CP 데이터베이스 (46M)
- **참고**: 데이터 변경 불가, 모든 프로젝트에서 공유

### **datasets/** - 중앙 데이터셋 저장소

- **용도**: 모든 실험의 입력 데이터
- **내용**:
  - `input_data/datasets/cp_positive_pairs.tsv`: CP 양성 쌍
  - `noncp_homolog_pairs.tsv`: 비-CP 호몰로그
  - `scope40_domains.tsv`: 도메인 메타데이터
- **설계**:
  - 중복 제거 (데이터 일관성 유지)
  - 모든 프로젝트가 동일 입력 사용 가능
  - 버전 관리 용이

### **experiments/** - 연구 프로젝트들

각 실험은 독립적인 프로젝트로, 다음 구조를 따름:
```
experiments/XX_프로젝트명/
├── README.md               # 프로젝트 설명 & 목표
├── scripts/                # 프로젝트 특화 스크립트
├── results/                # 분석 결과 (CSV, 그래프 등)
├── logs/                   # 실행 로그
└── notebooks/ (선택사항)   # Jupyter 분석
```

---

## 📊 **각 실험 요약**

### **01_encoder_training** 🔧 학습 모델
- **목표**: 3Di 인코더 학습 (8f, 9f, 10f)
- **상태**: ✅ 완료
- **주요 산출물**: encoder.pt, decoder.pt (encoders_and_tools/로 이동됨)

### **02_8f_10f_ssw_comparison** 📈 모델 비교
- **목표**: 8f vs 10f 모델의 SSW 정렬 성능 비교
- **상태**: ✅ 완료
- **주요 산출물**: 
  - ssw_comparison_8f_vs_10f.csv (267MB 대규모 비교)
  - 커버리지 플롯 (95% 커버리지 비교)

### **03_x1_x2_cp_detection** 🔍 CP 감지
- **목표**: X1/X2 정렬로 CP 검출 성능 평가
- **상태**: ✅ 완료
- **구조**:
  - `check_cp/`: CP 양성 쌍 분석 (451개)
  - `check_noncp/`: 비-CP 호몰로그 분석 (1,593개)
- **주요 결과**: noncp_8f_10f_comparison.csv

### **04_pdz_sanity_check** ✓ 검증 실험
- **목표**: PDZ 도메인 CP 감지 검증
- **상태**: ✅ 완료
- **실험**:
  1. Untrimmed SSW (원본 3Di 시퀀스)
  2. TM-align 구조 정렬 (기준선)
  3. Trimmed SSW (정렬 후 재비교)
- **주요 발견**:
  - 8f > 9f > 10f 순서로 CP 감지 우수
  - 8f의 AS 점수: 332 (다른 모델의 2-2.4배)

---

## 🔄 **데이터 흐름**

```
datasets/ (중앙)
    ↓
experiments/01 (모델 학습) → encoders_and_tools/
    ↓
experiments/02 (8f vs 10f)
    ↓
experiments/03 (X1/X2 CP 감지)
    ↓
experiments/04 (PDZ 검증)
```

---

## ⚠️ **중요 규칙**

### ✅ **유지해야 할 것**
- `core/`: 모델 & 도구 & 기본 데이터
- `datasets/`: 중앙 데이터셋 (모든 프로젝트가 참조)
- `experiments/XX_*/results/`: 최종 결과물
- `scripts/`: 공용 파이프라인 스크립트

### ❌ **삭제 가능**
- `tmp/`: 임시 작업 파일 (재현성 있으면 언제든 재생성)
- `experiments/XX_*/logs/`: 실행 로그
- 개별 실험의 중간 생성물

### 🔒 **공유 금지**
- `core/scope_databases/`: 데이터 변경 불가
- `core/encoders_and_tools/`: 모델 버전 잠금

---

## 📝 **활용 가이드**

### 새로운 실험 추가
```bash
# 1. 새 실험 디렉토리 생성
mkdir -p experiments/05_새프로젝트명
cd experiments/05_새프로젝트명

# 2. README 작성
cat > README.md << EOF
# 프로젝트명

## 목표
...

## 사용 데이터
- input_data/datasets/cp_positive_pairs.tsv
- core/encoders_and_tools/training_3di_gpu_*/

## 결과 위치
results/
EOF

# 3. 스크립트 작성
mkdir scripts results
```

### 기존 데이터 활용
```python
import pandas as pd
from pathlib import Path

# 중앙 데이터셋 로드
cp_pairs = pd.read_csv("../../input_data/datasets/cp_positive_pairs.tsv", sep='\t')

# 모델 경로
model_dir = Path("../../core/encoders_and_tools/training_3di_gpu_8f")
encoder = torch.load(model_dir / "encoder.pt")
```

---

## 📦 **로컬 셋업 체크리스트**

- [ ] `core/encoders_and_tools/` 모델 확인
- [ ] `core/scope_databases/` 데이터 확인
- [ ] `datasets/` 중앙 데이터셋 확인
- [ ] `scripts/pairwise_3di_pipeline.py` 실행 가능 확인
- [ ] 각 `experiments/XX_*/README.md` 읽기

---

## 🚀 **향후 계획**

- [ ] 추가 도메인 패밀리 분석 (protease, kinase 등)
- [ ] 8f 모델 특성 분석 (왜 최고 성능?)
- [ ] 9f 모델 재학습 (데이터셋 확대)
- [ ] 자동화 파이프라인 구축 (신규 도메인 자동 분석)
- [ ] 웹 대시보드 (결과 시각화)

---

**Last Updated**: 2025-12-15  
**Version**: 1.0  
**Maintainer**: @jugipalace
