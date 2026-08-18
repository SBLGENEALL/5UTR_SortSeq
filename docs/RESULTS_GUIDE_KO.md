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
| `expected_bin_score` | `6×p_bin1 + ... + 1×p_bin6` | 6에 가까울수록 high mCherry |
| `high15_probability` | `p_bin1 + p_bin2` | 해당 UTR의 gated cell이 top 15%에 있을 추정 확률 |
| `high15_enrichment` | high15 probability / 0.15 | 1=pool 평균, 1보다 크면 high 쪽 |
| `most_enriched_bin` | `p_ib / w_b`가 최대인 bin | pool 대비 가장 농축된 위치 |
| `dominant_cell_mass_bin` | `p_ib`가 최대인 bin | 실제 추정 cell mass가 가장 큰 위치 |
| `gate_entry_probability_capped` | whole unsorted 대비 target-gate representation 근사 | 보조 QC이며 주 발현 score가 아님 |
| `expression_tier` | score 기반 percentile tier | replicate 1개일 때 exact rank보다 권장 |

## 좋은 high 후보

- `pass_coverage=TRUE`
- `expected_bin_score` 상위 tier
- `high15_enrichment`가 충분히 큼
- `most_enriched_bin`이 bin1 또는 bin2
- bin1→bin6 확률이 한쪽으로 비교적 매끄럽게 이동
- 이상한 gate depletion이나 한-bin-only PCR jackpot이 없음
- 1-mismatch와 2-mismatch rescue 설정에서 tier가 안정적

## Original/reference 대비 우위

`REFERENCE_VARIANT_ID=auto`는 variant ID가 `original` 또는 `orginal`인 control을
대소문자와 관계없이 자동 탐지합니다. 다음 값은 `reference_comparison.tsv`와 주 결과표에
기록됩니다.

| 열 | 의미 |
|---|---|
| `delta_score_vs_reference` | UTR score − original score; 양수이면 high 쪽으로 이동 |
| `delta_high15_probability_vs_reference` | UTR의 bin1+2 probability − original 값 |
| `high15_fold_vs_reference` | UTR high15 probability / original high15 probability |
| `score_above_reference` | coverage 통과 및 score가 original보다 높음 |
| `score_and_high15_above_reference` | score와 high15 probability가 모두 original보다 높음 |

Score는 1–6 사이의 ordinal weighted average이므로 fold가 아니라 차이로 해석합니다.
High15 probability는 비율이므로 original 대비 fold로 비교할 수 있습니다. 그림에서는
original을 빨간 별, 빨간 테두리 또는 빨간 점선으로 표시합니다. Top50 밖에 있더라도
heatmap에는 비교 행으로 추가합니다.

기존 rescue/LibraryQC를 다시 계산하지 않고 reference 비교와 그림만 갱신하려면:

```bash
bash run_pipeline.sh reanalyze
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

## Strict coverage와 low-count 과대평가 방지

`pass_coverage`의 기본값(`unsorted >= 50`, `total six bins >= 100`)은 탐색 결과를
최대한 보존하기 위한 최소 기준입니다. 최종 hit에는 다음 strict 기준을 사용합니다.

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

## 결론 문구

biological replicate가 하나라면 다음 수준이 타당합니다.

> 이 UTR은 mCherry+/GFP- gate 내에서 높은 mCherry bin 쪽으로 이동했으며, cell-fraction-corrected score와 top-15% enrichment에서 상위 tier에 속해 후속 개별 construct 검증 후보로 선정하였다.

`translation rate를 증가시켰다`보다 `steady-state mCherry fluorescence가 높은 후보`라고 표현하세요. 정확한 1–2,000등, 통계적 FDR, 번역률의 직접 인과 결론은 independent biological replicate와 개별 construct 검증 전에는 피합니다.
