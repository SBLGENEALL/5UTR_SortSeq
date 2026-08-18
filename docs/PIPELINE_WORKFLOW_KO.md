# 5UTR Sort-seq 5단계 분석 흐름

이 문서는 시간이 지난 뒤에도 이번 분석이 **무엇을 위해, 어떤 순서로,
어떤 결과를 확인하며** 진행됐는지 재현할 수 있도록 만든 실행 기록입니다.

## 분석 환경과 원칙

- 분석 서버는 외부 인터넷이 차단된 **오프라인 Linux**입니다.
- GitHub ZIP은 인터넷 가능한 PC에서 받아 서버로 옮깁니다.
- 서버에서는 기존 pDNA QC Python/R 환경과 기존 `NGS_LibraryQC`를 재사용합니다.
- 서버에서 `curl`, `wget`, `git clone/pull`, 온라인 `pip`/R 설치를 하지 않습니다.
- 사내 raw/result 파일을 외부 서비스나 채팅에 업로드하지 않습니다. 문제 확인과 추가
  계산은 서버 내부 명령 및 로컬 HTML/TSV/PNG/PDF 산출물로 수행합니다.
- 원본 `raw_data`는 수정하지 않습니다.
- `bash run_pipeline.sh full`보다 아래 5단계를 하나씩 실행하고 결과를 확인하는
  방식을 권장합니다.

## 이번 실험의 입력과 방향

- sample: `bin1–bin6 + whole unsorted`, 총 7개
- `bin1`: mCherry가 가장 높은 구간
- `bin6`: mCherry-positive gate 안에서 가장 낮은 구간
- bin 내부 비율: `0.05, 0.10, 0.15, 0.20, 0.30, 0.20`
- whole unsorted의 `population_fraction`: 빈칸
- 전체 시작 세포의 50%를 unsorted로 따로 떼어낸 것은 bin 내부 비율에 넣지 않음
- target gate 전체 비율: 약 `0.90` (`OVERALL_GATE_FRACTION=0.90`)

## 분석 시작 전 한 번만 하는 설정

코드 폴더로 이동하고 기존 pDNA QC 환경을 활성화합니다. 실제 환경 활성화 명령은
서버에 이미 만들어 둔 환경에 맞게 사용합니다.

```bash
cd /data/user/MCET03/03_NGS/02_5UTR_sorting/5UTR_SortSeq-main
```

Run SampleSheet에서 실제 i7/i5를 자동으로 가져와 private config를 만듭니다.

```bash
python scripts/init_project.py \
  --raw-dir /data/user/MCET03/03_NGS/02_5UTR_sorting/raw_data \
  --config-dir /data/user/MCET03/03_NGS/02_5UTR_sorting/config \
  --run-sample-sheet /data/user/MCET03/03_NGS/02_5UTR_sorting/raw_data/Analysis/1/Data/260812_sample_sheet.csv
```

생성 파일:

```text
config/detected_fastq_samples.csv
config/SampleSheet.csv
config/sample_map.csv
```

실제 index가 채워졌는지 확인합니다.

```bash
if grep -q 'REPLACE_I7\|REPLACE_I5' \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/config/SampleSheet.csv
then
  echo 'ERROR: 실제 index가 추출되지 않음'
else
  echo 'OK: 실제 index 자동 추출 완료'
fi
```

실행 설정 경로와 CPU 설정을 한 번 지정합니다.

```bash
export SORTSEQ_PROJECT_CONFIG=/data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1
```

`sortseq.env`에서는 `WORKERS=128`, 기존 pDNA 분석에 실제로 사용한
`LIBRARYQC_PY`와 `LIBRARYQC_CONFIG`, raw/config/results 경로를 확인합니다.

현재 서버 성능과 실제 적용될 계산 설정도 확인합니다.

```bash
bash run_pipeline.sh resources
```

이번 서버에서 기대하는 주요 값은 다음과 같습니다.

| 항목 | 값/정책 |
|---|---|
| logical CPUs | 256 |
| physical cores | 128 |
| sockets / NUMA nodes | 2 / 4 |
| memory | 약 503 GiB |
| LibraryQC | 128 workers, batch size 10,000 |
| Rescue output | `pigz` 또는 내장 parallel Python gzip |
| NUMA | `numactl` 사용 가능할 때 memory interleave |

실제 탐지 결과는 `results/resource_profile.tsv`에도 기록됩니다.

## 전체 흐름 한눈에 보기

| 단계 | 명령 | 질문 | 핵심 산출물 |
|---|---|---|---|
| 1/5 | `preflight` | FASTQ와 dual index 설정이 맞는가? | `rescue_inspection.json` |
| 2/5 | `rescue` | Undetermined 중 어느 read를 안전하게 복구할 수 있는가? | rescued FASTQ, `rescue_manifest.json` |
| 3/5 | `libraryqc` | 각 sample에 어떤 UTR이 몇 read 있는가? | `variant_count_matrix.csv` |
| 4/5 | `analyze` | bin 크기와 NGS depth 보정 후 어떤 UTR이 high인가? | `utr_results_*.tsv` |
| 5/5 | `plot` | QC와 후보를 어떻게 시각적으로 검토할 것인가? | PNG, PDF |

## 실행 중 진행상황 확인

새 터미널에서 기존 분석 환경을 활성화하고 다음을 실행합니다.

```bash
cd /data/user/MCET03/03_NGS/02_5UTR_sorting/5UTR_SortSeq-main
export SORTSEQ_PROJECT_CONFIG=/data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env
bash run_pipeline.sh status
```

이 명령은 분석을 건드리지 않고 다음을 보여줍니다.

- 분석 프로세스가 `RUNNING`인지 `IDLE`인지
- 현재 또는 가장 최근 단계와 시작/완료/실패 상태
- rescue 진행률, 처리 read 수, 처리속도, 현재 rescue율과 ETA
- 1/5–5/5 checkpoint 파일 생성 여부

각 단계는 시작·완료 시각과 소요시간을 터미널에 표시하고 다음 파일에 남깁니다.

```text
results/pipeline_status.tsv
results/pipeline_history.tsv
```

Rescue는 기본 30초마다 다음과 비슷하게 갱신됩니다.

```text
[rescue] 42.50% | chunk 1/1 | processed 12,500,000 | rescued 10,000,000 (80.00%) | 25,000 read pairs/s | elapsed 8m 20s | ETA 11m 16s
```

진행률은 R1 gzip 파일에서 읽은 **압축 바이트 비율의 근사값**입니다. gzip read-ahead
때문에 아주 작게 흔들릴 수 있지만 대용량 FASTQ의 전체 진행상황과 ETA를 보는 데
적합합니다. 표시 간격은 private `sortseq.env`에서 바꿀 수 있습니다.

```bash
PROGRESS_INTERVAL_SECONDS=30
PROGRESS_CHECK_READS=100000
```

---

## 1/5. Preflight: 실제 분석 전 안전 점검

### 목적

- bin1–6, unsorted, Undetermined FASTQ 탐지
- FASTQ sample prefix와 run SampleSheet `Sample_ID`의 1:1 일치 확인
- FASTQ header의 `i7+i5` 존재 확인
- i5를 as-given으로 볼지 reverse complement로 볼지 자동 선택
- 1–2 mismatch rescue에 대한 expected dual-index 간 거리 확인

### 실행

```bash
bash run_pipeline.sh preflight 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/preflight.log
```

### 확인할 결과

```text
results/index_rescue/rescue_inspection.json
```

다음 조건을 확인합니다.

```text
missing_expected_samples_in_assigned_fastqs = []
unexpected_assigned_fastq_samples = []
undetermined_chunks >= 1
orientation_qc.headers_scanned_with_dual_index > 0
selected_i5_orientation = as_given 또는 reverse_complement
```

`detected_samples`에는 Undetermined가 있고 `expected_samples`에는 없는 것이
정상입니다. Undetermined는 8번째 생물학적 sample이 아니라 7개 expected index로
배정되지 못한 read 집합입니다.

---

## 2/5. Rescue: Undetermined dual-index 복구

### 목적

Undetermined FASTQ header의 i7+i5를 7개 expected dual-index pair와 비교합니다.

- exact, 1 mismatch, 2 mismatch 범위
- 가장 가까운 expected pair가 유일할 때만 해당 sample로 복구
- 동률, index 누락, 거리 초과, R1/R2 불일치는 residual로 보존
- 원본 FASTQ는 수정하지 않음

### 왜 큰 gzip 하나를 128개로 바로 나누지 않는가

현재 `undetermined_chunks=1`이므로 파일 단위 병렬 처리는 할 수 없습니다. 일반 gzip
스트림은 압축 바이트를 임의 위치에서 잘라 독립적으로 풀 수 없고, FASTQ도 반드시
4-line record와 R1/R2 pair 경계를 유지해야 합니다. 먼저 전부 해제해 수십 개의 임시
파일로 나눈 뒤 병렬 처리할 수는 있지만, 직렬 해제와 대용량 임시 파일 쓰기·읽기가
추가되어 보통 이 분석에는 불리합니다.

따라서 rescue는 입력을 한 번만 스트리밍하면서 다음을 병렬·최적화합니다.

- rescue 가능한 dual-index 조합을 미리 계산해 read당 7회 거리 계산을 O(1) lookup으로 변경
- sample별 R1/R2 출력 gzip을 독립 프로세스로 동시에 압축(최대 16 output streams)
- `pigz`가 서버에 이미 있으면 output당 기본 4 thread 사용
- `pigz`가 없어도 번들된 Python compressor process를 사용하므로 추가 설치 불필요

기본값은 private `sortseq.env`에 항목이 없어도 자동 적용됩니다.

```bash
RESCUE_COMPRESSION_BACKEND=auto
RESCUE_PIGZ_THREADS_PER_FILE=4
```

입력 gzip 해제 자체는 하나의 스트림이라는 제한이 남습니다. 다만 기존처럼 index 비교와
모든 결과 재압축까지 한 코어에서 처리하지 않기 때문에 rescue 전체 속도가 개선됩니다.
`rescue_manifest.json`의 `performance`에서 선택된 backend와 lookup hit 수를 확인할 수
있습니다.

### 실행

```bash
bash run_pipeline.sh rescue 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/rescue.log
```

실행 터미널에는 약 30초마다 진행률과 ETA가 출력되며, 다른 터미널에서는
`bash run_pipeline.sh status`로 같은 상태를 확인할 수 있습니다.

### 결과 구조

```text
results/index_rescue/
├── fastq/                         # 원본 symlink + rescued L900 chunk + residual
├── rescue_manifest.json
├── rescue_summary.csv
├── rescue_mismatch_distribution.csv
└── top_undetermined_barcodes.csv
```

Confidently rescued read는 bin1–6 또는 unsorted의 추가 FASTQ chunk에 포함됩니다.
`Undetermined_residual`은 다음 Library QC에서 별도 QC/count는 하지만 어느 bin인지
모르므로 최종 phenotype score에는 포함하지 않습니다.

### 완료 후 확인

```bash
cat /data/user/MCET03/03_NGS/02_5UTR_sorting/results/index_rescue/rescue_summary.csv
```

다음을 봅니다.

- 전체 Undetermined 중 rescue된 비율
- 특정 sample 하나로 지나치게 몰리지 않았는지
- 2-mismatch rescue가 0/1 mismatch보다 비정상적으로 많지 않은지
- ambiguous 또는 too-distant residual이 얼마나 남았는지

---

## 3/5. Library QC: FASTQ QC와 UTR counting

### 목적

기존 pDNA 분석과 같은 reference, primer/anchor, exact/near-match 규칙으로 모든
sample의 UTR raw count를 계산합니다. 이 단계가 `WORKERS=128`을 사용하는 주된
병렬 계산 단계입니다. 기본 batch size는 10,000이며, 503 GiB RAM에서 process 간
전달 횟수를 줄이도록 설정했습니다. 각 worker 내부 BLAS/OpenMP thread는 1개로 제한해
128 workers가 다시 중첩 thread를 생성하지 않게 합니다.

### 실행

```bash
bash run_pipeline.sh libraryqc 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/libraryqc.log
```

### 주요 결과

```text
results/library_qc/report.html
results/library_qc/RESULTS_TO_SHARE.txt
results/library_qc/combined/all_sample_metrics.csv
results/library_qc/combined/variant_count_matrix.csv
```

### 완료 후 확인

- 7개 phenotype sample이 모두 matrix column에 있는지
- rescued chunk가 원래 sample count에 합쳐졌는지
- sample별 assigned reference reads, Q30, anchor recovery, exact/near mapping
- UTR dropout과 count 불균형이 과도하지 않은지
- `Undetermined_residual`은 별도 QC column으로만 존재하는지

---

## 4/5. Analyze: bin 보정과 UTR phenotype score

### 목적

Raw read count를 직접 비교하지 않고 두 가지 차이를 보정합니다.

1. 각 NGS sample의 총 read depth 차이
2. 6개 sorter bin의 population 크기 차이

UTR `i`, bin `b`에 대해:

```text
within-bin frequency   f_ib = count_ib / total_assigned_reads_b
population mass        a_ib = bin_fraction_b × f_ib
bin probability        p_ib = a_ib / sum_b(a_ib)
```

### 실행

```bash
bash run_pipeline.sh analyze 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/analyze.log
```

### 주요 결과

```text
results/sortseq/utr_results_easy.tsv
results/sortseq/utr_results_full.tsv
results/sortseq/high_candidates.tsv
results/sortseq/strict_coverage_results.tsv
results/sortseq/high_confidence_candidates.tsv
results/sortseq/easy_report.html
```

### 결과를 보는 순서

1. `pass_coverage=TRUE`
2. `expected_bin_score`: 6에 가까울수록 높은 mCherry
3. `high15_probability`: bin1+bin2에 있을 추정 확률
4. `high15_enrichment`: 1이면 pool 평균, 1보다 크면 high 쪽 농축
5. `top15_vs_unsorted_log2_enrichment`: bin1+2의 whole-unsorted 대비 농축도;
   고발현 hit-selection의 주 랭킹
5. `most_enriched_bin`과 bin1–6 probability 분포
6. `expression_tier`: 정확한 개별 등수보다 top 1/5/10% tier 중심
7. whole unsorted 대비 target-gate representation은 보조 QC로 사용

최종 hit 목록은 permissive `pass_coverage`가 아니라 다음 strict filter를 추가로
통과한 UTR를 사용합니다.

```text
unsorted >= max(1,000, 전체 중앙값의 10%)
total six bins >= max(5,000, 전체 중앙값의 10%)
bin1 + bin2 raw support >= 200
3개 이상의 bin에서 검출
single-bin jackpot suspect 제외
```

Whole unsorted는 별도 기준 sample이므로 `population_fraction`을 비워 둡니다.
주 fluorescence score는 6개 bin 내부 분포에서 계산합니다.

### 계산을 한 단계씩 직접 감사하기

LibraryQC의 `variant_count_matrix.csv`를 사용해 전체 UTR의 중간 계산값을 CSV로
내보낼 수 있습니다.

```bash
bash run_pipeline.sh scoring-steps
```

기존 rescue/LibraryQC를 유지하고 v0.1.9 top15 랭킹만 다시 계산할 때는 환경을
나누어 실행합니다.

```bash
# Python 환경
bash run_pipeline.sh reanalyze-analysis

# R 환경
bash run_pipeline.sh plot
```

Python 단계는 기존 `results/sortseq`만 `archive/`로 이동합니다. `raw_data`,
`results/index_rescue`, `results/library_qc`는 이동하거나 다시 계산하지 않습니다.
계산식과 후보 flag는 [TOP15_ENRICHMENT_KO.md](TOP15_ENRICHMENT_KO.md)를 참조하세요.

| 파일 | 계산 단계 | 의미 |
|---|---|---|
| `01_raw_counts.csv` | Raw count | variant count matrix의 원래 count |
| `02_depth_normalized_frequency.csv` | Depth normalization | `c_ib / N_b`; 각 bin 내 UTR 비율 |
| `03_bin_size_corrected_mass.csv` | Bin-size correction | `w_b × c_ib/N_b` |
| `04_within_utr_bin_probability.csv` | Within-UTR normalization | UTR별 bin1–6 합을 1로 정규화 |
| `05_score_contributions_and_final_score.csv` | Weighted score | probability에 6–1점을 곱해 합산 |
| `06_all_steps_combined_audit.csv` | 전체 결합 | 한 UTR의 모든 중간값을 한 행에서 추적 |

`probability`는 NGS read 확률이 아니라 다음 조건부 세포분율 추정치입니다.

```text
P(bin b | UTR i, target gate)
= w_b × (count_ib / total_reads_b)
  / sum_k[w_k × (count_ik / total_reads_k)]
```

Unsorted는 raw/depth-normalized QC에는 포함되지만 6-bin fluorescence score에는
포함되지 않습니다.

---

## 5/5. Plot: R QC figure와 PDF

### 목적

Index rescue, sample depth, UTR score, high-bin enrichment, top 후보 분포,
unsorted 대비 target gate를 한 번에 검토할 그림을 만듭니다.

### 실행

```bash
bash run_pipeline.sh plot 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/plot.log
```

### 결과

```text
results/sortseq/figures/01_undetermined_rescue.png
results/sortseq/figures/02_sample_assigned_reads.png
results/sortseq/figures/03_score_vs_high15_enrichment.png
results/sortseq/figures/04_top_utr_bin_heatmap.png
results/sortseq/figures/05_expression_tiers.png
results/sortseq/figures/06_unsorted_vs_target_gate.png
results/sortseq/figures/sortseq_qc_figures.pdf
results/sortseq/figures/strict/08_strict_score_vs_high15_enrichment.png
results/sortseq/figures/strict/09_strict_top_utr_bin_probability_heatmap.png
results/sortseq/figures/strict/10_strict_top_utr_bin_enrichment_heatmap.png
results/sortseq/figures/strict/11_strict_unsorted_vs_target_gate.png
results/sortseq/figures/strict/12_strict_reference_score_position.png
results/sortseq/figures/strict/strict_candidate_figures.pdf
```

09번은 절대 cell fraction이고, 10번은 `log2(probability / bin fraction)`입니다.
따라서 bin3처럼 원래 더 큰 bin이 09번에서 밝게 보이더라도, 10번에서는 실제 상대
농축이 bin1 또는 bin2인지 바로 확인할 수 있습니다.

## 중단 후 재개하는 방법

완료된 단계를 다시 실행하지 않고 다음 단계부터 이어갑니다.

```text
preflight 완료 → rescue 실행
rescue 완료    → libraryqc 실행
libraryqc 완료 → analyze 실행
analyze 완료   → plot 실행
```

예를 들어 rescue까지 완료됐다면 다음부터 시작합니다.

```bash
bash run_pipeline.sh libraryqc
bash run_pipeline.sh analyze
bash run_pipeline.sh plot
```

이미 rescue 결과가 있는데 `full`을 실행하면 기존 결과 보호 기능 때문에 중단될 수
있습니다. `full --replace`는 기존 결과를 archive로 옮기고 처음부터 다시 계산할 때만
사용합니다.

## 최종 해석의 한계

- biological replicate가 하나이면 `estimated_rank`는 탐색적 순위입니다.
- 정확한 1–2,000등이나 통계적 FDR보다 top tier 후보 선정이 타당합니다.
- 결과는 `steady-state mCherry fluorescence가 높은 후보`로 표현합니다.
- `translation rate를 직접 증가시켰다`는 결론에는 개별 construct와 독립 replicate,
  필요 시 mRNA/translation 분해 측정이 추가로 필요합니다.
