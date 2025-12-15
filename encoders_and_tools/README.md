# Encoders and Tools

## 구성

```
encoders_and_tools/
├── training_3di_gpu_8f/      # 8f VQ-VAE 인코더
├── training_3di_gpu_10f/     # 10f Foldseek 관련 파일 (참고용)
├── foldseek_10f/             # Foldseek 10f 메트릭 및 구성
└── ssw/                       # SSW (Smith-Waterman) 정렬 도구
```

## 주요 파일

### training_3di_gpu_8f/
- `tmp/encoder.pt`: 8f VQ-VAE 인코더 모델
- `tmp/states.txt`: 8f 알파벳 상태 파일
- `encode_pdbs.py`: PDB → 8f 3Di 변환 스크립트

### foldseek_10f/
- `s_10f.mat`: 10f 대체 매트릭스

### ssw/
- `tmp/ssw/src/ssw_test`: SSW 실행 바이너리
- `s_8f.mat`: 8f 대체 매트릭스

## 사용법

### 8f 인코딩
```bash
cd training_3di_gpu_8f
echo "domain1.pdb" | python encode_pdbs.py tmp/encoder.pt tmp/states.txt --virt 270 0 2
```

### SSW 정렬
```bash
ssw/tmp/ssw/src/ssw_test -a matrix.mat -p -c query.fasta target.fasta
```

## 참고

- 이 디렉토리의 도구들은 모든 분석에서 공유됨
- 상대 경로로 접근할 때는 필요에 따라 `../encoders_and_tools/` 사용
