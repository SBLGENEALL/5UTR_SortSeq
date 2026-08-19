# 결과 해석 가이드

## 핵심 계산

UTR `i`, bin `b`의 count를 `c_ib`, 해당 sample의 전체 assigned read를 `N_b`, sorter population fraction을 `w_b`라고 하면:

```text
within-bin frequency      f_ib = c_ib / N_b
bin-corrected mass        a_ib = w_b × f_ib
UTR별 bin probability     p_ib = a_ib / sum_b(a_ib)
```

따라서 raw count가 큰 bin을 고르는 것이 아니라, sequencing depth와 실제 bin 크기를 모두 보정한 `p_ib`를 사용합니다.

## `utr_results_easy.tsv`에서 보는 순서

| 열 | 의미 | 권장 해석 |
|---|---|---|
| `pass_coverage` | unsorted와 6-bin count cutoff 통과 | 먼저 `TRUE`만 봄 |
| `high15_final_rank` | technical-bootstrap robust High15 순위 | **cloning primary; 1이 우선** |
| `high15_probability` | `p_bin1 + p_bin2` | 해당 UTR의 gated cell이 top 15%에 있을 추정 확률; primary point estimate |
| `expected_bin_score` | `6×p_bin1 + ... + 1×p_bin6` | 전체 분포 high 이동의 supporting metric |
| `high15_enrichment` | high15 probability / 0.15 | 1=pool 평균, 1보다 크면 high 쪽 |
| `high15_robust_rank_score` | bootstrap `log2(High15/reference High15)`의 10th percentile | 낮은-read 우연을 보수적으로 낮춘 ranking 값 |
| `high15_bootstrap_probability_above_comparator` | resample 중 High15가 reference보다 높은 비율 | NGS read-sampling 안정성; biological probability 아님 |
| `candidate_tier` | tier1/tier2/tier3/not_candidate | cloning 우선순위와 review 경고 |
| `top15_vs_unsorted_log2_enrichment` | combined bin1+2 frequency / unsorted frequency의 log2 | gate representation 포함 보조/탐색 지표 |
| `bin1_vs_unsorted_enrichment` | bin1 frequency / unsorted frequency | extreme-high 구간 농축도 |
| `bin2_vs_unsorted_enrichment` | bin2 frequency / unsorted frequency | 두 번째 high 구간 농축도 |
| `most_enriched_bin` | `p_ib / w_b`가 최대인 bin | pool 대비 가장 농축된 위치 |
| `dominant_cell_mass_bin` | `p_ib`가 최대인 bin | 실제 추정 cell mass가 가장 큰 위치 |
| `gate_entry_probability_capped` | whole unsorted 대비 target-gate representation 근사 | 보조 QC이며 주 발현 score가 아님 |
| `expression_tier` | score 기반 percentile tier | replicate 1개일 때 exact rank보다 권장 |

## 좋은 high 후보

- `pass_coverage=TRUE`
- `high15_final_rank`가 높고 `high15_probability`가 original보다 큼
- `high15_technical_stability_pass=TRUE`
- `candidate_tier`가 tier1 또는 tier2
- tier1이면 `expected_bin_score`도 original보다 큼
- `most_enriched_bin`이 bin1 또는 bin2
- bin1→bin6 확률이 한쪽으로 비교적 매끄럽게 이동
- 이상한 gate depletion이나 한-bin-only PCR jackpot이 없음
- 1-mismatch와 2-mismatch rescue 설정에서 tier가 안정적

고발현 hit 선정은 `high15_final_rank`를 주 랭킹으로 사용합니다. 전체 분포와
intermediate-high phenotype은 expected score로 보조 확인하고, unsorted enrichment는
gate representation QC/탐색값으로 분리합니다. 계산과 후보 기준은
[TOP15_ENRICHMENT_KO.md](TOP15_ENRICHMENT_KO.md)에 정리되어 있습니다.

## Original/reference 대비 우위

`REFERENCE_VARIANT_ID=auto`는 variant ID가 `original` 또는 `orginal`인 control을
대소문자와 관계없이 자동 탐지합니다. 다음 값은 `reference_comparison.tsv`와 주 결과표에
기록됩니다.

| 열 | 의미 |
|---|---|
| `delta_score_vs_reference` | UTR score − original score; 양수이면 high 쪽으로 이동 |
| `delta_high15_probability_vs_reference` | UTR의 bin1+2 probability − original 값 |
| `high15_fold_vs_reference` | UTR high15 probability / original high15 probability |
| `high15_log2_fold_vs_reference` | 위 비율의 log2; fluorescence fold가 아님 |
| `score_above_reference` | coverage 통과 및 score가 original보다 높음 |
| `score_and_high15_above_reference` | score와 high15 probability가 모두 original보다 높음 |

Score는 1–6 사이의 ordinal weighted average이므로 fold가 아니라 차이로 해석합니다.
High15 probability는 비율이므로 original 대비 fold로 비교할 수 있습니다. 그림에서는
original을 빨간 별, 빨간 테두리 또는 빨간 점선으로 표시합니다. Top50 밖에 있더라도
heatmap에는 비교 행으로 추가합니다.

기존 rescue/LibraryQC를 다시 계산하지 않고 Python 분석과 R 그림을 나누어 갱신하려면:

```bash
bash run_pipeline.sh reanalyze-analysis
# R 환경으로 전환한 뒤
bash run_pipeline.sh plot
```

## unsorted의 역할

현재 unsorted는 whole population이고 6 bins는 mCherry+/GFP- gate 내부입니다. 따라서 unsorted는 주 fluorescence score의 분모가 아닙니다. plasmid/unsorted representation과 target-gate 진입의 보조 QC로 사용합니다.

`gate_entry_probability_raw > 1`은 실제 확률이 100%를 넘었다는 뜻이 아니라 composition/PCR/sampling 차이를 확인하라는 QC 신호입니다.

## `bin_probability`의 정확한 의미

`bin1_probability`–`bin6_probability`는 raw-read 비율이 아닙니다. 다음 두 보정을 한 뒤
해당 UTR 내부에서 합이 1이 되도록 정규화한 값입니다.

1. sequencing depth: `count_ib / total_assigned_UTR_reads_b`
2. sorter bin size: 위 값에 `population_fraction_b`를 곱함

따라서 `bin1_probability=0.20`은 해당 UTR를 가진 target-gate 세포 중 약 20%가
bin1에 존재한다고 추정한다는 뜻입니다. 직접 관찰한 단일세포 확률이나 transfection
확률은 아닙니다. 계산 중간값은 `bash run_pipeline.sh scoring-steps`로 생성되는
`results/sortseq/scoring_steps/` CSV에서 확인합니다.

## Read support와 low-count 과대평가 방지

`pass_coverage`의 기본값(`unsorted >= 50`, `total six bins >= 100`)은 탐색 결과를
최대한 보존하기 위한 최소 기준입니다. 기존 score 중심 strict 표는 다음 강한-support
QC를 유지하지만, 새 High15 cloning rank의 hard filter는 아닙니다.

| 열 | 기본 기준 |
|---|---:|
| `unsorted_count` | `max(1,000, 중앙값의 10%)` 이상 |
| `total_6bin_count` | `max(5,000, 중앙값의 10%)` 이상 |
| `high_bin_raw_count` | 200 이상 |
| `detected_in_n_bins` | 3 이상 |

`strict_coverage_pass=TRUE`이면서 기존 high-score/high15 조건을 통과하고,
`single_bin_jackpot_suspect=FALSE`인 UTR만
`high_confidence_candidate_flag=TRUE`가 됩니다. 결과는 다음 파일에 기록됩니다.

```text
results/sortseq/strict_coverage_results.tsv
results/sortseq/high_confidence_candidates.tsv
```

R 그림은 `results/sortseq/figures/strict/`에 생성됩니다. Probability heatmap은 절대
세포분율이므로 큰 bin이 밝게 보일 수 있습니다. Enrichment heatmap은
`log2(probability / population fraction)`을 사용하므로 bin1–6의 서로 다른 크기를
제거하고 어느 bin에 상대적으로 농축됐는지 보여줍니다.

새 cloning list는 `total_6bin_count >= 200`, `unsorted_count >= 50`,
`high_bin_raw_count >= 20`을 eligibility로 사용한 뒤 technical bootstrap과 jackpot
flag로 안정성을 구분합니다. 핵심 파일은 `top_candidates_for_cloning.csv`입니다.

## UTR A형 양극화 분포

bin1·2와 bin5·6가 동시에 높고 bin3·4가 낮은 UTR는 평균적인 고발현 UTR와 구분합니다.
이 분포는 일부 cell이 high일 가능성을 보여주지만, pooled 6-bin 데이터만으로 두 개의
single-cell peak를 증명하지는 못합니다. 또한 낮은 raw count에서는 우연히 middle bin이
비어 양극화처럼 보일 수 있습니다.

```bash
bash run_pipeline.sh bimodality-qc
```

결과의 `clear_bimodal_flag`는 total six-bin count 200 이상, 양쪽 tail raw count 각각
20 이상, high/low tail과 middle-valley 조건을 모두 통과한 exploratory flag입니다.
`bimodality_by_read_count.csv`와 `bimodality_count_sensitivity.csv`에서 cutoff를 높일수록
후보가 급감하는지 확인합니다. 자세한 설명은
[BIMODALITY_QC_KO.md](BIMODALITY_QC_KO.md)를 보세요.

## 두 scoring endpoint 비교

`expected_bin_score`와 `high15_probability`를
`bash run_pipeline.sh compare-metrics`로 total six-bin count가 200보다 큰 동일 UTR
집합에서 Spearman rho와 Top50 overlap으로 진단합니다. High15가 primary이므로
consensus를 필수 조건으로 삼지 않으며, score-only와 top15-only의 분포 차이를 24번
heatmap에서 확인합니다. 자세한 실행과 판단은
[METRIC_COMPARISON_KO.md](METRIC_COMPARISON_KO.md)를 보세요.

## 결론 문구

biological replicate가 하나라면 다음 수준이 타당합니다.

> 이 UTR은 단일 pooled sort에서 mCherry+/GFP- gate 내부의 High15 probability가 original보다 높고 기술적 read-bootstrap에서 안정적으로 유지되어, 개별 clone 검증 후보로 우선 선정하였다. Six-bin weighted score는 전체 분포의 high 이동을 보조적으로 지지하였다.

`translation rate를 증가시켰다`보다 `steady-state mCherry fluorescence가 높은 후보`라고 표현하세요. 정확한 1–2,000등, 통계적 FDR, 번역률의 직접 인과 결론은 independent biological replicate와 개별 construct 검증 전에는 피합니다.
