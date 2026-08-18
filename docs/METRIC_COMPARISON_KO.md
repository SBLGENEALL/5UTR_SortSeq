# Expected score와 bin1+bin2 enrichment 비교

## 비교 목적

두 지표는 서로 틀린 계산이 아니라 서로 다른 phenotype을 요약합니다.

| 지표 | 질문 |
|---|---|
| `expected_bin_score` | UTR 세포의 전체 6-bin 분포가 평균적으로 high 쪽인가? |
| `top15_vs_unsorted_log2_enrichment` | whole unsorted 대비 bin1+bin2 extreme-high 구간에 농축됐는가? |

이번 비교에서는 UTR마다 6개 bin raw count 합이 **200 이하이면 제외**합니다. 코드에서는
`METRIC_COMPARE_MIN_TOTAL_COUNT=201`로 설정합니다.

Top15 enrichment는 unsorted를 분모로 사용하므로, 최종 비교 universe에는 추가로
다음 최소 read support를 적용합니다.

```text
total six-bin count >= 201
unsorted count >= 50
bin1+bin2 raw count >= 20
두 지표 모두 finite
```

`total-count-only`와 위 `robust-read-support` 조건의 상관성을 모두 계산하므로 unsorted
저coverage가 결과를 흔드는지도 확인할 수 있습니다.

## 실행

Python 환경:

```bash
cd /data/user/MCET03/03_NGS/02_5UTR_sorting/5UTR_SortSeq-main
export SORTSEQ_PROJECT_CONFIG=/data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env

bash run_pipeline.sh compare-metrics 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/compare_metrics.log
```

기존 rescue, LibraryQC, analyze, bimodality 결과는 수정하지 않습니다.

R 환경:

```bash
bash run_pipeline.sh plot 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/plot_metric_comparison.log
```

## Python 결과

```text
results/sortseq/metric_comparison/metric_comparison_summary.csv
results/sortseq/metric_comparison/metric_comparison_summary.json
results/sortseq/metric_comparison/metric_comparison_all.csv
results/sortseq/metric_comparison/metric_comparison_supported.csv
results/sortseq/metric_comparison/metric_comparison_topn_overlap.csv
results/sortseq/metric_comparison/top_candidates_consensus.csv
results/sortseq/metric_comparison/top_candidates_score_only.csv
results/sortseq/metric_comparison/top_candidates_top15_only.csv
```

- `consensus`: 두 방식 모두 Top 50
- `score_only`: expected score에서만 Top 50
- `top15_only`: bin1+2 enrichment에서만 Top 50

## R 결과

```text
results/sortseq/figures/metric_comparison/21_expected_score_vs_top15_enrichment.png
results/sortseq/figures/metric_comparison/22_score_rank_vs_top15_rank.png
results/sortseq/figures/metric_comparison/23_top_candidate_overlap.png
results/sortseq/figures/metric_comparison/24_discordant_candidate_bin_heatmap.png
results/sortseq/figures/metric_comparison/metric_comparison_figures.pdf
```

21번은 실제 metric 값, 22번은 순위 관계, 23번은 Top 50 중복, 24번은 한 방식에서만
선정된 후보의 6-bin 분포를 보여줍니다.

## 가장 먼저 볼 값

터미널에는 다음 값이 바로 출력됩니다.

```text
Shared comparison universe
Spearman rho
Top 50 overlap
Agreement class
```

Spearman rho는 정확한 직선 관계보다 두 순위의 전반적인 일치도를 평가합니다.

| 결과 | 탐색적 해석 |
|---|---|
| rho ≥ 0.7, Top50 overlap ≥ 60% | 강한 일치; 상위 후보 대부분 공유 |
| rho 0.4–0.7 | 부분 일치; discordant 후보 분포 확인 필요 |
| rho < 0.4 | 서로 다른 phenotype을 선택하는 경향 |

위 숫자는 통계적 또는 생물학적 cutoff가 아니라 의사결정을 쉽게 하기 위한 기술적
가이드입니다. p-value보다 rho 크기와 Top50 overlap을 우선 봅니다.

## 어떤 지표를 최종 사용해야 하는가

실험 목표가 **일부 세포라도 extreme-high bin1+2에 들어가는 UTR**이라면
`top15_vs_unsorted_log2_enrichment`를 주 지표로 사용합니다. 다만 top15-only 후보가
bin5·6에도 많으면 heterogeneous high-tail 후보이므로 균일한 고발현으로 부르지 않습니다.

목표가 **세포 전체가 전반적으로 high 쪽에 위치하는 안정적인 UTR**이라면
`expected_bin_score`를 주 지표로 사용합니다. 다만 score-only 후보가 bin3 중심이면
extreme-high보다는 intermediate-high 후보입니다.

현재 목적처럼 “mCherry가 높은 5′UTR”를 처음 선별할 때는 다음 순서를 권장합니다.

1. `top_candidates_consensus.csv`를 우선 검증
2. original보다 두 지표가 모두 높은 후보를 우선
3. top15 enrichment를 extreme-high 주 endpoint로 사용
4. expected score와 bin5+6 probability를 안정성/heterogeneity 보조 지표로 사용
5. score-only와 top15-only는 24번 heatmap으로 별도 검토

## 분석 방식의 문헌 연결

- Matreyek et al. (2018): 여러 fluorescence bin의 weighted average로 전체 phenotype
  위치를 요약하는 방식. <https://pmc.ncbi.nlm.nih.gov/articles/PMC5980760/>
- Cao et al. (2021): 상위 fluorescence bin의 unsorted 대비 enrichment로 high-expression
  sequence를 선별하는 방식. <https://www.nature.com/articles/s41467-021-24436-7>
- Peterman and Levine (2016): threshold enrichment와 mean fluorescence가
  cell-to-cell variability에 따라 다른 결과를 낼 수 있음을 설명.
  <https://link.springer.com/article/10.1186/s12864-016-2533-5>

