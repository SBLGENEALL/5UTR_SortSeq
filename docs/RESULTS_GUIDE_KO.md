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

## unsorted의 역할

현재 unsorted는 whole population이고 6 bins는 mCherry+/GFP- gate 내부입니다. 따라서 unsorted는 주 fluorescence score의 분모가 아닙니다. plasmid/unsorted representation과 target-gate 진입의 보조 QC로 사용합니다.

`gate_entry_probability_raw > 1`은 실제 확률이 100%를 넘었다는 뜻이 아니라 composition/PCR/sampling 차이를 확인하라는 QC 신호입니다.

## 결론 문구

biological replicate가 하나라면 다음 수준이 타당합니다.

> 이 UTR은 mCherry+/GFP- gate 내에서 높은 mCherry bin 쪽으로 이동했으며, cell-fraction-corrected score와 top-15% enrichment에서 상위 tier에 속해 후속 개별 construct 검증 후보로 선정하였다.

`translation rate를 증가시켰다`보다 `steady-state mCherry fluorescence가 높은 후보`라고 표현하세요. 정확한 1–2,000등, 통계적 FDR, 번역률의 직접 인과 결론은 independent biological replicate와 개별 construct 검증 전에는 피합니다.
