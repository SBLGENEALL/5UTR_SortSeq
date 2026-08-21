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
# analyze가 Step 1-8 export, Step 6-8 비교, profile QC까지 자동 실행
bash run_pipeline.sh bimodality-qc
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
results/sortseq/strict_coverage_results.tsv
results/sortseq/high_confidence_candidates.tsv
results/sortseq/reference_comparison.tsv
results/sortseq/high15_primary_ranking.csv
results/sortseq/high15_candidates.csv
results/sortseq/top_candidates_for_cloning.csv
results/sortseq/top200_candidates_for_cloning.csv
results/sortseq/top15_enrichment_ranking.csv
results/sortseq/top15_candidates.csv
results/sortseq/top15_priority_candidates.csv
results/sortseq/all_utr_profiles/all_utr_profile_assignments.csv
results/sortseq/all_utr_profiles/all_utr_profile_cluster_summary.csv
results/sortseq/all_utr_profiles/all_utr_normalized_enrichment_profiles.csv
results/sortseq/all_utr_profiles/top200_normalized_enrichment_profiles.csv
results/sortseq/bimodality_qc/bimodality_summary.csv
results/sortseq/bimodality_qc/clear_bimodal_candidates.csv
results/sortseq/bimodality_qc/utra_like_strong_polarization.csv
results/sortseq/metric_comparison/metric_comparison_summary.csv
results/sortseq/metric_comparison/step6_step7_step8_correlations.csv
results/sortseq/metric_comparison/step6_step7_step8_topn_overlap.csv
results/sortseq/metric_comparison/top_candidates_consensus.csv
results/sortseq/figures/sortseq_qc_figures.pdf
results/sortseq/figures/strict/strict_candidate_figures.pdf
results/sortseq/figures/top15/top15_candidate_figures.pdf
results/sortseq/figures/all_utr_profiles/all_utr_profile_figures.pdf
```

`original` 또는 `orginal` control UTR는 자동 탐지됩니다. Reference 대비 High15와
보조 score 차이가 표에 추가되고, scatter/heatmap/score-distribution/unsorted-gate 그림에는
빨간 별 또는 선으로 표시됩니다. Rescue와 LibraryQC를 유지하고 Python 분석만
다시 만들려면 다음을 실행합니다.

```bash
bash run_pipeline.sh reanalyze-analysis

# 별도 R 환경에서
bash run_pipeline.sh plot
```

기존 `results/sortseq`은 `archive/`로 이동되며 raw FASTQ, rescue, LibraryQC 결과는
수정되지 않습니다.

## High15 primary ranking과 cloning 후보

v0.2.0부터 cloning 후보의 1차 지표는 target gate 내부의
`high15_probability = P(bin1 또는 bin2 | UTR, mCherry+/GFP-)`입니다.
`expected_bin_score`는 전체 6-bin 이동을 확인하는 보조 지표이고,
`top15_vs_unsorted_log2_enrichment`는 gate representation이 포함된 탐색/QC 지표입니다.

기본 phenotype eligibility는 `total six bins >= 200`, `bin1+bin2 >= 20`입니다.
Unsorted는 conditional High15의 분모나 hard filter가 아니라 representation 및 gate-entry
QC로 별도 보고합니다.
각 bin의 variant count vector를 1,000회 기술적 read-bootstrap하여
`log2(High15/original High15)`의 10th percentile로 robust rank를 만듭니다. 이 bootstrap은
PCR 또는 biological uncertainty를 추정하지 않습니다. bin1과 bin2의 개별 unsorted
enrichment는 hard filter로 사용하지 않습니다.

기본 cloning/plot shortlist는 Top 200입니다. `top200_candidates_for_cloning.csv`와
normalized-enrichment profile 그림에는 `original`의 순위와 분포가 반드시 표시됩니다.

최종 `top_candidates_for_cloning.csv`에는 안정적인 tier1 clean high-shift와 tier2
high-tail 후보가 High15 순으로 들어갑니다. 자세한 식과 열 해석은
[High15 중심 고발현 랭킹](docs/TOP15_ENRICHMENT_KO.md)을 보세요.

## UTR A형 high+low-tail 분포 QC

bin1·2와 bin5·6에 동시에 많고 bin3·4가 비어 보이는 UTR는 평균 score만으로 분류하지
않습니다. `bash run_pipeline.sh bimodality-qc`는 양쪽 tail probability, middle valley,
양쪽 tail의 raw-read support를 함께 검사하고, total read 구간별 후보 비율과
read count–polarization 상관을 계산합니다. 이 결과는 “bimodality-like shape” QC이며
두 biological state의 증명은 아닙니다.

Python 환경에서 QC를 실행한 뒤 R 환경에서 기존 `plot` 명령을 실행하면
`results/sortseq/figures/bimodality/`에 17–20번 그림과 PDF가 생성됩니다. 자세한 기준과
해석은 [UTR A형 양극화 분포 QC](docs/BIMODALITY_QC_KO.md)를 보세요.

## Step 6·7·8 직접 비교

`bash run_pipeline.sh compare-metrics`는 동일한 read-supported UTR 집합에서 Step 6
conditional High15, Step 7 expected score, Step 8 `p/w` relative enrichment와
`q=f/sum(f)=(p/w)/sum(p/w)` normalized-enrichment profile을 계산합니다.
Step 6–7의 Spearman 상관성과 Top 20/100/200 overlap을 보고하고, Step 6–8은 rho=1,
overlap=100%인지 자동 검산합니다. R `plot`을 실행하면
`figures/metric_comparison/`에 21–29번 그림이 생성됩니다. 자세한 해석은
[두 scoring 방식 직접 비교](docs/METRIC_COMPARISON_KO.md)를 보세요.

## 전체 UTR의 6-bin 분포 보기

Top 후보뿐 아니라 read-support 기준을 통과한 모든 UTR의 분포를 함께 보려면 Python
환경과 R 환경에서 차례로 실행합니다.

```bash
# Python 환경
bash run_pipeline.sh profile-qc

# R 환경
bash run_pipeline.sh plot
```

`profile-qc`는 모든 UTR에 대해 `q=f/sum(f)=(p/w)/sum(p/w)`를 계산하여
`all_utr_normalized_enrichment_profiles.csv`에 저장합니다. Cluster와 전체-library 그림은
**bin1–bin6 raw count 합이 200 이상**인 UTR에 적용하며, `bin1+2` 또는 unsorted cutoff는
사용하지 않습니다. 별도로 High15 Top 200과 `original`을 heatmap 및 25개씩 나눈
profile page로 표시합니다. 정확한 cluster와 전체 UTR 값은
`all_utr_profile_assignments.csv` 및 `all_utr_normalized_enrichment_profiles.csv`에서
확인합니다. 자세한 설명은
[전체 UTR 6-bin 분포 시각화](docs/ALL_UTR_PROFILE_QC_KO.md)를 보세요.

## 계산 전체를 Step 1–8 CSV로 확인

`variant_count_matrix.csv`에서 바로 시작해 각 중간 계산을 확인하려면 Python 환경에서
다음을 실행합니다.

```bash
bash run_pipeline.sh scoring-steps
```

이 명령은 기존 분석 결과를 삭제하지 않고 raw count부터 profile enrichment까지 모든
중간 계산 CSV를 추가합니다. `analyze`를 새로 실행한 경우에는 자동으로 생성됩니다.

```text
results/sortseq/scoring_steps/00_sample_parameters.csv
results/sortseq/scoring_steps/01_raw_counts.csv
results/sortseq/scoring_steps/02_depth_normalized_frequency.csv
results/sortseq/scoring_steps/03_bin_size_corrected_mass.csv
results/sortseq/scoring_steps/04_within_utr_bin_probability.csv
results/sortseq/scoring_steps/05_score_contributions_and_final_score.csv
results/sortseq/scoring_steps/06_all_steps_combined_audit.csv
results/sortseq/scoring_steps/07_calculation_checks.csv
results/sortseq/scoring_steps/08_high15_primary_ranking.csv
results/sortseq/scoring_steps/09_top15_unsorted_enrichment_secondary.csv
results/sortseq/scoring_steps/step01_raw_counts.csv
results/sortseq/scoring_steps/step02_depth_normalized_frequency.csv
results/sortseq/scoring_steps/step03_bin_size_corrected_mass.csv
results/sortseq/scoring_steps/step04_corrected_mass_total.csv
results/sortseq/scoring_steps/step05_within_utr_probability.csv
results/sortseq/scoring_steps/step06_high15_probability.csv
results/sortseq/scoring_steps/step07_expected_score.csv
results/sortseq/scoring_steps/step08_relative_enrichment_profile.csv
results/sortseq/scoring_steps/step09_all_steps_audit.csv
```

여기서 `binN_probability`는 sequencing read가 그 bin에 들어갈 확률이 아니라,
depth와 bin size를 보정한 뒤 추정한
`P(bin N | 해당 UTR, mCherry+/GFP- gate)`입니다. 한 UTR의 bin1–6 probability 합은
1이며, 이 probability에 `6,5,4,3,2,1`을 곱한 합이 expected bin score입니다.
Step 8의 `binN_relative_enrichment=p_ib/w_b`는 probability가 아니라 bin 기본 크기 대비
농축도이며, `log2(p_ib/w_b)`가 0이면 중립입니다.

## 기존 strict score QC와 보조 그림

기본 `pass_coverage`는 전체 탐색표를 보존하기 위한 느슨한 기준입니다. 기존 score
중심 분석의 강한 read-support QC에는 다음 strict filter가 유지됩니다.

```text
unsorted count >= max(1,000, 전체 UTR 중앙값의 10%)
total six-bin count >= max(5,000, 전체 UTR 중앙값의 10%)
bin1 + bin2 raw count >= 200
검출 bin 수 >= 3
한두 bin에 probability 85% 이상 집중된 경우 jackpot 의심 표시
```

`bash run_pipeline.sh analyze`는 `strict_coverage_results.tsv`와
`high_confidence_candidates.tsv`를 만듭니다. R 환경에서 `bash run_pipeline.sh plot`을
실행하면 기존 그림을 유지하면서 `figures/strict/` 아래에 strict 후보용 그림과 PDF를
추가합니다. 사용자가 서버에서 이미 `high_confidence_candidates.tsv`를 만든 경우에는
analyze를 다시 실행하지 않고 plot만 실행해도 됩니다. 파일이 없거나 비어 있어도
v0.1.6부터는 `utr_results_full.tsv`에서 같은 strict 기준을 R이 자동 재계산합니다.

관측된 count 중앙값이 고정 strict floor보다 훨씬 낮아 strict high 후보가 0개이면,
v0.1.7의 R 그림은 `unsorted >= 100`, `total six bins >= 500`, `bin1+2 >= 50`,
`detected bins >= 3`을 만족하는 UTR를 **read-supported exploratory** 후보로 표시합니다.
`unsorted >= 200`, `total >= 1,000`, `bin1+2 >= 100`인 strong-support UTR는 별도
색으로 표시하며, strong-support high 후보가 0개라는 사실도 그대로 보존합니다.

이 strict score 표는 v0.2.0의 High15 cloning primary rank를 대체하지 않습니다.
실제 cloning 우선순위는 `top_candidates_for_cloning.csv`와 `candidate_tier`를 봅니다.

특히 두 heatmap을 구분해서 봅니다.

- `09_strict_top_utr_bin_probability_heatmap.png`: UTR별 절대 bin probability
- `10_strict_top_utr_bin_enrichment_heatmap.png`: `log2(probability / bin fraction)`;
  bin 크기를 기준으로 상대 농축된 곳은 빨강, depletion된 곳은 파랑

전체 흐름은 [docs/PIPELINE_WORKFLOW_KO.md](docs/PIPELINE_WORKFLOW_KO.md),
서버 설정은 [docs/SERVER_GUIDE_KO.md](docs/SERVER_GUIDE_KO.md), 결과 해석은
[docs/RESULTS_GUIDE_KO.md](docs/RESULTS_GUIDE_KO.md), Undetermined 처리 원칙은
[docs/UNDETERMINED_RESCUE_KO.md](docs/UNDETERMINED_RESCUE_KO.md), 양극화 분포는
[docs/BIMODALITY_QC_KO.md](docs/BIMODALITY_QC_KO.md), scoring 선택은
[docs/METRIC_COMPARISON_KO.md](docs/METRIC_COMPARISON_KO.md)를 보세요.

## 자체 테스트

```bash
python -m unittest discover -s tests -v
bash -n run_pipeline.sh
```

현재 biological replicate가 하나라면 `estimated_rank`는 탐색적 순위입니다. 정확한 1–2,000등이나 FDR보다 top 1/5/10% tier와 개별 construct 재검증을 사용하세요.
