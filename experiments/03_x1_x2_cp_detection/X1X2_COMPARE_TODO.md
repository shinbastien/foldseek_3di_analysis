# x1_x2_compare 디렉토리 정리 상태

## 📋 현재 상황 (2025-12-15)

### 현재 Root에 존재하는 디렉토리들
```
x1_x2_compare/           (648K) - 주 결과 파일들
x1_x2_compare_test/      (640K) - 테스트용 파일들  
x1_x2_compare_noncp_random/ (336K) - 비교용 데이터
```

### ⚠️ TODO: Experiment 03 완료 후 정리

**이 파일들을 `experiments/03_x1_x2_cp_detection/results/`로 이동할 예정입니다.**

현재는 root에 유지되고 있으며, experiment 03의 데이터 정리가 완료된 후에 다음을 진행해야 합니다:

1. **확인 필요사항**
   - [ ] x1_x2_compare_test/: 실제로 테스트용인지, 아니면 중요한 데이터인지 확인
   - [ ] 각 디렉토리의 내용 및 용도 재확인

2. **이동 계획**
   - [ ] x1_x2_compare/ → experiments/03_x1_x2_cp_detection/results/main_results/
   - [ ] x1_x2_compare_test/ → experiments/03_x1_x2_cp_detection/results/test_results/ (또는 삭제)
   - [ ] x1_x2_compare_noncp_random/ → experiments/03_x1_x2_cp_detection/results/noncp_comparison/

3. **완료 후 cleanup**
   - root 레벨에서 제거

---

**작성일**: 2025-12-15  
**우선순위**: Experiment 03 정리 완료 후 처리
