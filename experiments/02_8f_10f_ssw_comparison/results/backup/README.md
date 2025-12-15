# 8f vs 10f SSW Score Comparison

## 개요
8f (커스텀 VQ-VAE) 과 10f (Foldseek) 3Di 알파벳을 직접 비교하는 분석입니다.

## 데이터 구조

```
8f_10f_ssw_comparison/
├── results/                          # 최종 결과
│   ├── comparison_results.csv       # 메인 결과
│   ├── metrics_summary.csv          # 통계 요약
│   └── plots/                        # 시각화
│       └── comparison_plots.png
├── tmp/                              # 중간 파일
├── logs/                             # 실행 로그
│   └── processing.log
└── scripts/                          # 실행 스크립트
    ├── compare_8f_10f_ssw.py        # 메인 스크립트
    └── compare_8f_10f_ssw_backup.py # (백업)
```

## 실행 방법

```bash
cd 8f_10f_ssw_comparison

python scripts/compare_8f_10f_ssw.py \
  [필요한 인자들]
```

(자세한 사용법은 스크립트 내 주석 참고)

## 결과

- `results/comparison_results.csv`: 8f vs 10f 점수 비교
- `results/metrics_summary.csv`: 통계 요약
- `results/plots/comparison_plots.png`: 시각화

## 참고

- 관련 인코더/도구는 `encoders_and_tools/` 폴더에 위치
