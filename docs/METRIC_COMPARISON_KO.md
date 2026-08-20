# Step 6·7·8 비교와 후보 순위 해석

## 세 단계가 나타내는 것

| 단계 | 계산 | 의미 | 역할 |
|---|---|---|---|
| Step 6 | `H_i = p_i1 + p_i2` | target gate 안에서 해당 UTR 세포가 상위 15%에 있을 추정 비율 | **primary cloning rank** |
| Step 7 | `S_i = 6p_i1+5p_i2+...+p_i6` | 전체 6-bin 분포의 평균적인 high/low 위치 | supporting endpoint |
| Step 8 | `R_ib=p_ib/w_b`, `L_ib=log2(R_ib)` | 각 bin의 원래 크기보다 해당 UTR가 얼마나 농축되었는지 | profile 시각화 |

Step 8의 high-tail 요약값은 다음과 같습니다.

```text
Step8_High15_i = High15_i / (w1+w2) = High15_i / 0.15
```

따라서 Step 6과 Step 8 high-tail 값은 숫자의 단위만 다르고 **순위는 완전히
같아야 합니다**. Spearman rho=1, Top20/50/100 overlap=100%, rank mismatch=0이
정상입니다. 다르게 나오면 생물학적 차이가 아니라 계산 또는 필터 적용 오류입니다.

실제로 정보를 주는 비교는 Step 6과 Step 7입니다.

- 둘 다 높은 UTR: 상위 tail이 높고 전체 분포도 high 쪽으로 이동한 clean high-shift
- Step 6만 높은 UTR: bin1·2 tail은 높지만 middle/low mass도 있는 tail-heavy 후보
- Step 7만 높은 UTR: 전체적으로 위로 이동했지만 bin1·2에 특별히 집중되지는 않은 후보

두 지표의 Top-N 교집합만을 hard filter로 쓰지는 않습니다. 최종 cloning 우선순위는
Step 6으로 고정하고 Step 7과 Step 8 profile을 보조 근거로 검토합니다.

## 공통 비교 대상

기본값은 다음과 같습니다.

```text
six-bin total raw count >= 200
bin1+bin2 raw count >= 20
Step 6·7·8 값이 모두 finite
```

Unsorted는 Step 6·7·8 계산의 분모가 아니므로 기본 cutoff는 0입니다. Unsorted count는
library representation 및 gate-entry QC에 별도로 남습니다. 필요하면
`METRIC_COMPARE_MIN_UNSORTED_COUNT`를 양수로 설정해 추가 QC universe를 만들 수 있습니다.

## 실행 순서

`analyze`는 Python 환경에서 Step 1–8 export, Step 6–8 비교, 전체 profile clustering까지
자동 실행합니다.

```bash
# Python 환경
bash run_pipeline.sh analyze 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/analyze.log

# R 환경
bash run_pipeline.sh plot 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/plot.log
```

이미 `analyze`가 끝난 결과에서 비교만 다시 계산하려면 다음 명령을 사용합니다.

```bash
bash run_pipeline.sh compare-metrics 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/compare_metrics.log
```

이 명령들은 FASTQ, rescue 및 LibraryQC 결과를 삭제하지 않습니다.

## Python 결과

```text
results/sortseq/metric_comparison/step6_step7_step8_metrics.csv
results/sortseq/metric_comparison/step6_step7_step8_correlations.csv
results/sortseq/metric_comparison/step6_step7_step8_topn_overlap.csv
results/sortseq/metric_comparison/step6_step7_step8_summary.json
```

이전 버전과의 호환을 위해 다음 이름도 계속 생성됩니다.

```text
metric_comparison_all.csv
metric_comparison_supported.csv
metric_comparison_summary.csv
metric_comparison_topn_overlap.csv
top_candidates_consensus.csv
top_candidates_score_only.csv
top_candidates_top15_only.csv
```

`step6_step7_step8_topn_overlap.csv`의 `pair` 값:

- `step6_vs_step7`: primary와 supporting score의 실제 후보 겹침
- `step6_vs_step8`: 항상 100%여야 하는 계산 self-check
- `step7_vs_step8`: Step 7과 primary-equivalent Step 8의 겹침
- `step6_step7_step8_three_way`: 세 목록의 교집합; Step 6=8이므로 Step 6·7 교집합과 같음

## R 그림

```text
21_expected_score_vs_high15_probability.png
22_score_rank_vs_top15_rank.png
23_top_candidate_overlap.png
24_discordant_candidate_bin_heatmap.png
25_step6_step7_step8_correlation.png
26_step6_step7_step8_topn_overlap.png
27_top_high15_relative_enrichment_heatmap.png
28_group_median_relative_enrichment_profiles.png
metric_comparison_figures.pdf
```

27·28번은 x축을 `bin6 → bin1`로 배치합니다. 따라서 고발현 쪽으로 농축된 UTR는
오른쪽이 올라가거나 붉게 나타납니다. y축은 `log2(p/w)`입니다.

- `0`: 그 bin의 nominal size와 같은 중립 수준
- `+1`: nominal 대비 2배 농축
- `-1`: nominal 대비 절반으로 depletion

절대 probability `p`를 그리는 그림과 달리, 큰 bin3·5가 원래 크기 때문에 밝아지는
효과를 제거하고 profile의 방향을 보여줍니다.

## 결과 판단

1. `step6_step7_step8_summary.json`에서 `step8_identity_check_pass=true` 확인
2. Step 6–7 Spearman rho로 전체 rank 경향 확인
3. Top20/50/100 overlap으로 실제 상위 후보 겹침 확인
4. `top_candidates_consensus.csv`는 clean high-shift 보조 목록으로 사용
5. `top_candidates_top15_only.csv`도 폐기하지 말고 tail-heavy 후보로 개별 검증
6. 최종 primary 순위는 Step 6 High15로 유지

Biological replicate가 없으므로 correlation p-value나 candidate p-value/FDR을
생물학적 유의성으로 해석하지 않습니다. 이 비교는 단일 pooled sort 안에서 후보 선정
방식의 일관성과 profile을 점검하는 기술적·서술적 분석입니다.
