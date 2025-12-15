# Root 레벨 파일 정리 상태

## 📋 현재 상황 (2025-12-15)

### 이동된 파일들
- ✓ `pairwise_3di_pipeline.py` → `experiments/01_pdz_sanity_check/scripts/`
- ✓ `tmalign_3di_match_pipeline.py` → `experiments/01_pdz_sanity_check/scripts/`
- ✓ `run_pipeline_batch.sh` → `experiments/01_pdz_sanity_check/batches/`

### Root에 유지된 이유
Root 레벨의 파일들은 아직 이동되지 않은 상태입니다. 다음 사항을 고려해야 합니다:

1. **경로 참조 호환성**: 이 스크립트들이 다른 파일/패키지에서 import되거나 실행될 수 있음
2. **레거시 호환성**: 기존 실행 명령어나 스크립트가 이들을 root에서 찾을 수 있음
3. **Utility 함수 이동**: `utils/pdz_pipeline_utils.py`로 공통 함수 추출 완료

## ⚠️ 다음 단계 필요

Root 레벨의 파일들을 완전히 정리하려면:

1. **원본 파일 삭제 여부 결정**
   - Option A: Root 레벨 원본 삭제 (복사본만 experiment에 유지)
   - Option B: Root 레벨 유지 (백업/호환성 목적)

2. **Script 경로 업데이트**
   - 스크립트 내부의 import 경로 확인
   - `utils/pdz_pipeline_utils.py` import 추가
   - 상대 경로 조정 필요 가능성

## 📝 Action Items

다음 것들을 확인/실행해야 함:
- [ ] experiment 버전의 스크립트들이 `utils/pdz_pipeline_utils.py`를 올바르게 import하는지 확인
- [ ] Root 레벨 원본 파일들을 삭제할지 유지할지 결정
- [ ] 필요시 symbolic link 생성 고려

---
**작성일**: 2025-12-15
