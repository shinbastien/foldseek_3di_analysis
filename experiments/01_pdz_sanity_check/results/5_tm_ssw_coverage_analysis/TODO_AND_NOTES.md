# 5_tm_ssw_coverage_analysis - TODO & 상황 정리

## 📋 현재 상황

이 디렉토리는 PDZ 도메인 구조 비교 실험의 **5번 단계 (TM-align + SSW coverage analysis)** 결과를 저장합니다.

### 디렉토리 구조
```
5_tm_ssw_coverage_analysis/
├── tmp_pipeline_outputs/           # tmalign_3di_match_pipeline.py의 중간 결과 (2024-11-25 생성)
│   ├── 2hga_trim_vs_2vsv_PDZ_w{0,1,2}/
│   ├── 2vsv_PDZ_vs_2z9i_pdz_w{0,1,2}/
│   ├── 2z9i_pdz_vs_2hga_trim_w{0,1,2}/
│   └── Sap_C_vs_Sap_C_circular_permutation_w{0,1,2}/
└── [최종 결과 파일들 - 아직 정리 안됨]
```

## ⚠️ 주요 문제점

### 1. **재생성 가능성 미확인**
   - `tmp_pipeline_outputs/` 폴더들은 `tmalign_3di_match_pipeline.py`를 실행하면 재생성될 수 있음
   - 각 폴더 내 `logs/runner.out` 파일을 보면 실행 명령어가 기록되어 있음
   - **확인 필요**: 동일한 결과를 생성하는가? (변수/시드 때문에 다를 수 있음)

### 2. **결과 파일 분산**
   - TM_align 결과 (`.txt`): `tmp_pipeline_outputs/` 각 폴더 내 존재
   - 최종 요약 CSV: 아직 정리되지 않음
   - **정리 필요**: 최종 결과 파일을 `results/` 디렉토리로 수합

### 3. **공간 효율성**
   - `tmp_pipeline_outputs/`: ~상당한 크기 (정확한 크기 미측정)
   - **고려사항**: 
     - 필요한 결과만 추출 후 tmp 폴더 삭제할지
     - 아니면 완전한 재현성을 위해 모든 중간 결과 유지할지

## 📝 TODO

### 우선순위 1 (즉시)
- [ ] `tmp_pipeline_outputs/` 각 폴더에서 최종 결과 파일 (`.csv`, `.txt`) 추출
- [ ] 결과 파일을 `5_tm_ssw_coverage_analysis/` 직하위로 이동
- [ ] 최종 요약 파일 생성

### 우선순위 2 (재생성 테스트 후)
- [ ] `tmalign_3di_match_pipeline.py` 재실행하여 결과 재현성 확인
- [ ] 재현 가능하면 `tmp_pipeline_outputs/` 삭제 권장
- [ ] 재현 불가능하면 현재 위치 유지

### 우선순위 3 (문서화)
- [ ] 각 폴더의 `runner.out`에서 실행 명령어 추출
- [ ] 재실행 스크립트 작성 (필요시)

## 🔗 관련 스크립트
- `../scripts/tmalign_3di_match_pipeline.py`: 이 단계의 메인 스크립트
- `../scripts/pairwise_3di_pipeline.py`: 4번 단계의 스크립트
- `../../utils/pdz_pipeline_utils.py`: 공유 utility 함수

## 📊 재생성 명령어 예시 (runner.out에서 추출)

각 `tmp_pipeline_outputs/*/logs/runner.out` 파일을 참고하여:
1. foldseek createdb
2. pdb_to_3di.py (8f, 9f, 10f variants)
3. TMalign 실행
4. SSW 실행
5. 결과 파일 생성

---

**작성일**: 2025-12-15  
**상태**: 진행 중 - 최종 결과 파일 정리 및 재생성성 테스트 대기
