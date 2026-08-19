# High15 probability와 expected score 비교

## 두 지표의 역할

| 지표 | 질문 | 현재 역할 |
|---|---|---|
| `high15_probability` | target gate 안의 해당 UTR 세포 중 bin1+2에 있을 비율은? | **primary cloning rank** |
| `expected_bin_score` | 해당 UTR의 전체 6-bin 분포가 평균적으로 얼마나 high 쪽인가? | supporting metric |

두 지표의 후보가 다르게 나오는 것은 계산 오류가 아닙니다. High15는 상위 tail을,
score는 전체 분포의 평균 위치를 요약합니다. 현재 실험 목표가 bin1+2의 mCherry-high
UTR 선택이므로 Top-N overlap을 필수 조건으로 사용하지 않습니다.

## 공통 비교 대상

```text
total six-bin count >= 201
unsorted count >= 50
bin1+bin2 raw count >= 20
두 지표 모두 finite
```

Unsorted cutoff는 primary High15 식의 분모로 쓰기 위한 것이 아니라, 동일한 기존
read-support universe와 gate-representation QC를 유지하기 위한 eligibility 조건입니다.

## 실행

Python 환경:

```bash
bash run_pipeline.sh compare-metrics 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/compare_metrics.log
```

R 환경:

```bash
bash run_pipeline.sh plot 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/plot_metric_comparison.log
```

기존 FASTQ, rescue, LibraryQC, Python 분석표는 이 명령으로 삭제되지 않습니다.

## 결과

```text
results/sortseq/metric_comparison/metric_comparison_summary.csv
results/sortseq/metric_comparison/metric_comparison_summary.json
results/sortseq/metric_comparison/metric_comparison_supported.csv
results/sortseq/metric_comparison/metric_comparison_topn_overlap.csv
results/sortseq/metric_comparison/top_candidates_consensus.csv
results/sortseq/metric_comparison/top_candidates_score_only.csv
results/sortseq/metric_comparison/top_candidates_top15_only.csv

results/sortseq/figures/metric_comparison/21_expected_score_vs_high15_probability.png
results/sortseq/figures/metric_comparison/22_score_rank_vs_top15_rank.png
results/sortseq/figures/metric_comparison/23_top_candidate_overlap.png
results/sortseq/figures/metric_comparison/24_discordant_candidate_bin_heatmap.png
```

- `consensus`: 두 방식 모두 Top N
- `score_only`: expected score에서만 Top N
- `top15_only`: conditional High15 probability에서만 Top N

Spearman rho와 Top-N overlap은 두 요약이 얼마나 비슷한지 보여주는 진단값입니다.
높거나 낮다고 해서 primary endpoint를 자동으로 바꾸지 않습니다.

## 후보 선택에 사용하는 방법

1. `high15_final_rank`로 전체 우선순위를 정합니다.
2. `expected_bin_score > original`이면 clean high-shift를 지지하여 tier1로 표시합니다.
3. score가 original 이하라도 High15가 안정적이면 tier2 high-tail 후보로 유지합니다.
4. 한 bin 몰빵, 낮은 bootstrap stability는 tier3 review로 내립니다.
5. `top_candidates_for_cloning.csv`에서 tier1·tier2를 우선 검증합니다.

따라서 과거의 “두 지표 Top50 교집합만 선택” 규칙은 사용하지 않습니다. 교집합만
고집하면 서로 다른 phenotype 요약을 동시에 극단화하여 유용한 high-tail 후보를
과도하게 잃을 수 있습니다.

Biological replicate가 없으므로 상관 p-value, candidate p-value 또는 FDR을 생물학적
유의성으로 해석하지 않습니다.
