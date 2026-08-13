# 5UTR_SortSeq

`whole unsorted + mCherry high→low 6 bins`로 구성된 7-sample 5′UTR Sort-seq 분석 파이프라인입니다.

이 저장소는 기존 [`NGS_LibraryQC`](https://github.com/SBLGENEALL/NGS_LibraryQC)를 대체하지 않습니다. 역할을 다음처럼 분리합니다.

| 단계 | 담당 |
|---|---|
| FASTQ/amplicon QC, reference UTR counting | `NGS_LibraryQC` |
| Undetermined dual-index rescue | `5UTR_SortSeq` |
| 6개 bin 크기·NGS depth 보정 | `5UTR_SortSeq` |
| UTR score/tier/QC/R plot | `5UTR_SortSeq` |

## 이번 실험에서 고정한 bin 방향

| sample | 의미 | nominal fraction | score |
|---|---|---:|---:|
| bin1 | 가장 높은 mCherry | 0.05 | 6 |
| bin2 | 다음 high | 0.10 | 5 |
| bin3 |  | 0.15 | 4 |
| bin4 |  | 0.20 | 3 |
| bin5 |  | 0.30 | 2 |
| bin6 | mCherry-positive gate 내 가장 낮음 | 0.20 | 1 |
| unsorted | gate 전 whole population | 해당 없음 | 주 score에 사용하지 않음 |

## 서버에서 가장 짧은 실행법

```bash
cd /data/user/MCET03/03_NGS/02_5UTR_sorting
unzip /옮겨놓은/5UTR_SortSeq.zip
cd 5UTR_SortSeq

python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
Rscript scripts/install_r_packages.R
```

FASTQ 이름에서 7개 sample prefix를 자동으로 읽어 설정 초안을 만듭니다.

```bash
python scripts/init_project.py \
  --raw-dir /data/user/MCET03/03_NGS/02_5UTR_sorting/raw_data \
  --config-dir /data/user/MCET03/03_NGS/02_5UTR_sorting/config
```

예를 들어 `UTR_bin1_A2_s1_R1_001.fastq.gz`는 `UTR_bin1_A2`로 인식됩니다. 생성된 `SampleSheet.csv`의 `REPLACE_I7`, `REPLACE_I5`를 실제 7개 index로 바꾸고 `sample_map.csv`를 확인합니다.

```bash
cp config/project.env.example \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env

# sortseq.env의 NGS_LibraryQC script/config 경로를 확인한 뒤 실행
export SORTSEQ_PROJECT_CONFIG=/data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env
bash run_pipeline.sh preflight
bash run_pipeline.sh full
```

`preflight`에서 index가 header에 있는지, i5 방향, sample 이름, index 간 최소 거리를 먼저 확인합니다. 기존 결과를 보존하며 다시 실행하려면 `bash run_pipeline.sh full --replace`를 사용합니다. 이전 결과는 `archive/`로 이동합니다.

## 먼저 볼 결과

```text
results/index_rescue/rescue_inspection.json
results/index_rescue/rescue_summary.csv
results/library_qc/report.html
results/library_qc/combined/variant_count_matrix.csv
results/sortseq/utr_results_easy.tsv
results/sortseq/utr_results_full.tsv
results/sortseq/high_candidates.tsv
results/sortseq/figures/sortseq_qc_figures.pdf
```

상세 실행법은 [docs/SERVER_GUIDE_KO.md](docs/SERVER_GUIDE_KO.md), 결과 해석은 [docs/RESULTS_GUIDE_KO.md](docs/RESULTS_GUIDE_KO.md), Undetermined 처리 원칙은 [docs/UNDETERMINED_RESCUE_KO.md](docs/UNDETERMINED_RESCUE_KO.md)를 보세요.

## 자체 테스트

```bash
python -m unittest discover -s tests -v
bash -n run_pipeline.sh
```

현재 biological replicate가 하나라면 `estimated_rank`는 탐색적 순위입니다. 정확한 1–2,000등이나 FDR보다 top 1/5/10% tier와 개별 construct 재검증을 사용하세요.
