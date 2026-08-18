# 5UTR_SortSeq

`whole unsorted + mCherry high→low 6 bins`로 구성된 7-sample 5′UTR Sort-seq 분석 파이프라인입니다.

> **권장 실행 방식:** 전체를 한 번에 실행하기보다
> [5단계 분석 흐름과 단계별 명령](docs/PIPELINE_WORKFLOW_KO.md)에 따라
> `preflight → rescue → libraryqc → analyze → plot`을 하나씩 실행하고 각 결과를
> 확인하세요. 이 문서는 각 단계의 목적, 완료 판정, 중단 후 재개 방법도 설명합니다.

> **서버 환경:** 분석 서버는 오프라인 Linux입니다. GitHub ZIP은 인터넷 가능한
> PC에서 받아 서버로 옮기고, 서버에서는 기존 pDNA QC Python/R 환경과
> `NGS_LibraryQC`를 재사용합니다.
>
> **데이터 보안:** 사내 FASTQ, count table, TSV/CSV 및 결과 파일은 외부로 업로드하거나
> 반출하지 않습니다. 코드만 GitHub ZIP으로 반입하고 모든 계산·QC·그림 생성은 사내
> 오프라인 Linux 서버 안에서 수행합니다.

실행 중 진행상황은 별도 터미널에서 다음 명령으로 확인할 수 있습니다.

```bash
export SORTSEQ_PROJECT_CONFIG=/data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env
bash run_pipeline.sh status
```

Rescue 단계는 30초마다 진행률, 처리 read 수, 속도, rescue율, 경과시간과 ETA를
출력하고 `results/index_rescue/rescue_progress.json`에도 기록합니다.

서버의 128 physical core/256 logical CPU/4 NUMA node/503 GiB 구성을 파이프라인이
어떻게 사용할지 먼저 확인할 수 있습니다.

```bash
bash run_pipeline.sh resources
```

기본 성능 설정은 `NGS_LibraryQC` 128 workers, batch size 10,000, NUMA memory
interleave입니다. Undetermined rescue는 index lookup을 미리 계산하고 각 sample의
R1/R2 gzip 출력을 별도 프로세스로 동시에 압축합니다. 서버에 `pigz`가 이미 있으면
자동 사용하고, 없어도 내장 Python 병렬 압축을 사용하므로 온라인 설치는 필요 없습니다.

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

## 서버에서 시작하기

```bash
cd /data/user/MCET03/03_NGS/02_5UTR_sorting
unzip /옮겨놓은/5UTR_SortSeq.zip
cd 5UTR_SortSeq-main
```

기존 pDNA QC에 사용한 Python/R 환경을 활성화합니다. 오프라인 서버에서 온라인
`pip`, R package 설치, `curl`, `wget`, `git clone/pull`을 실행하지 않습니다.

FASTQ 이름에서 7개 sample prefix를 읽고, 원래 NGS run SampleSheet의 실제
`Sample_ID/index/index2`를 검증하여 설정 파일을 만듭니다.

```bash
python scripts/init_project.py \
  --raw-dir /data/user/MCET03/03_NGS/02_5UTR_sorting/raw_data \
  --config-dir /data/user/MCET03/03_NGS/02_5UTR_sorting/config \
  --run-sample-sheet /data/user/MCET03/03_NGS/02_5UTR_sorting/raw_data/Analysis/1/Data/260812_sample_sheet.csv
```

`--run-sample-sheet`를 생략해도 `raw_data/Analysis` 아래에 sample-sheet CSV가
정확히 하나뿐이면 자동으로 찾습니다. 예를 들어
`UTR_bin1_A2_s1_R1_001.fastq.gz`는 `UTR_bin1_A2`로 인식되고, 원래
SampleSheet의 같은 `Sample_ID`에서 실제 i7/i5를 가져옵니다. 생성된
`SampleSheet.csv`에는 index가 이미 채워지므로 `sample_map.csv`의 bin 방향만
확인하면 됩니다.

```bash
cp config/project.env.example \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env

# sortseq.env의 NGS_LibraryQC script/config 경로를 확인한 뒤 단계별 실행
export SORTSEQ_PROJECT_CONFIG=/data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env
bash run_pipeline.sh preflight
bash run_pipeline.sh rescue
bash run_pipeline.sh libraryqc
bash run_pipeline.sh analyze
bash run_pipeline.sh plot
```

각 명령이 끝날 때마다 결과를 확인한 후 다음 단계로 넘어갑니다. 자세한 완료 판정은
[5단계 분석 흐름](docs/PIPELINE_WORKFLOW_KO.md)을 따르세요. 처음부터 완전히
재분석해야 할 때만 `bash run_pipeline.sh full --replace`를 사용합니다.

## 먼저 볼 결과

```text
results/index_rescue/rescue_inspection.json
results/index_rescue/rescue_summary.csv
results/resource_profile.tsv
results/library_qc/report.html
results/library_qc/combined/variant_count_matrix.csv
results/sortseq/utr_results_easy.tsv
results/sortseq/utr_results_full.tsv
results/sortseq/high_candidates.tsv
results/sortseq/reference_comparison.tsv
results/sortseq/figures/sortseq_qc_figures.pdf
```

`original` 또는 `orginal` control UTR는 자동 탐지됩니다. Reference 대비 score 차이와
high15 fold가 표에 추가되고, scatter/heatmap/score-distribution/unsorted-gate 그림에는
빨간 별 또는 선으로 표시됩니다. Rescue와 LibraryQC를 유지하고 reference-aware 결과만
다시 만들려면 새 코드에서 다음을 실행합니다.

```bash
bash run_pipeline.sh reanalyze
```

기존 `results/sortseq`은 `archive/`로 이동되며 raw FASTQ, rescue, LibraryQC 결과는
수정되지 않습니다.

전체 흐름은 [docs/PIPELINE_WORKFLOW_KO.md](docs/PIPELINE_WORKFLOW_KO.md),
서버 설정은 [docs/SERVER_GUIDE_KO.md](docs/SERVER_GUIDE_KO.md), 결과 해석은
[docs/RESULTS_GUIDE_KO.md](docs/RESULTS_GUIDE_KO.md), Undetermined 처리 원칙은
[docs/UNDETERMINED_RESCUE_KO.md](docs/UNDETERMINED_RESCUE_KO.md)를 보세요.

## 자체 테스트

```bash
python -m unittest discover -s tests -v
bash -n run_pipeline.sh
```

현재 biological replicate가 하나라면 `estimated_rank`는 탐색적 순위입니다. 정확한 1–2,000등이나 FDR보다 top 1/5/10% tier와 개별 construct 재검증을 사용하세요.
