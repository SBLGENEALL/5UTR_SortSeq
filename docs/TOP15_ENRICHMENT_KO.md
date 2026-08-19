# High15 중심 고발현 5′UTR 후보 랭킹

## 결론

이번 파이프라인의 cloning 후보 **1차 순위**는 다음 값입니다.

```text
High15 probability = P(bin1 또는 bin2 | UTR, mCherry+/GFP- target gate)
```

`expected_bin_score`는 전체 6-bin 분포가 high 방향으로 이동했는지 확인하는
**보조 지표**입니다. Whole unsorted를 분모로 한 enrichment는 target-gate
representation까지 섞이므로 primary rank가 아니라 별도 QC/탐색 지표로 남깁니다.

## 계산식

UTR `i`, bin `b`의 raw count를 `c_ib`, bin 전체 assigned read를 `N_b`, sorter
population fraction을 `w_b`라 정의합니다.

```text
1. Depth normalization       f_ib = c_ib / N_b
2. Bin-size correction       a_ib = w_b × f_ib
3. Within-UTR normalization  p_ib = a_ib / Σ_k a_ik
4. Primary endpoint          H_i  = p_i1 + p_i2
5. Supporting score          S_i  = 6p_i1 + 5p_i2 + 4p_i3 + 3p_i4 + 2p_i5 + p_i6
```

현재 `w=(0.05, 0.10, 0.15, 0.20, 0.30, 0.20)`이며, 따라서 pool-neutral
High15는 `0.15`입니다. `original/orginal`이 검출되면 pool 0.15 대신 reference의
실제 `H_original`과 비교합니다.

```text
high15_fold_vs_reference      = H_i / H_original
high15_log2_fold_vs_reference = log2(H_i / H_original)
```

이 비율은 fluorescence의 배수 증가가 아닙니다. “reference보다 top 15% 구간에
들어갈 추정 비율이 몇 배인가”라는 뜻입니다.

## Unsorted 지표를 분리하는 이유

```text
A_i = Σ_b w_b f_ib
G_i = A_i / u_i
E_i = top15_vs_unsorted_enrichment
    = G_i × H_i / 0.15
```

`E_i`에는 순수한 high-tail 위치 `H_i`뿐 아니라 whole unsorted 대비 target-gate
representation `G_i`도 포함됩니다. 이 실험에서 whole unsorted는 동일 시점·동일
세포군의 유용한 QC이지만, 목표가 target gate 내부의 mCherry-high UTR 선택이므로
`H_i`를 primary로 사용합니다. `G_i`와 `E_i`는 결과에서 삭제하지 않습니다.

## Read-support eligibility

기본값은 다음과 같습니다.

```text
unsorted raw count >= 50
six-bin total raw count >= 200
bin1 + bin2 raw count >= 20
```

이는 dropout과 극단적으로 희소한 UTR를 제외하는 실용적 eligibility filter이며,
통계적 유의성이나 정밀도를 보장하는 cutoff가 아닙니다.

## 기술적 bootstrap과 최종 순위

기본 1,000회 bootstrap에서 각 bin의 전체 variant count vector를 원래 read depth로
multinomial resampling하고, 위 1–4단계를 다시 계산합니다.

```text
high15_robust_rank_score
  = bootstrap log2(H_i / H_original)의 10th percentile

high15_final_rank
  = high15_robust_rank_score 내림차순 순위
```

함께 제공되는 값:

- `high15_bootstrap_probability_above_comparator`: resample 중 `H_i>H_original`인 비율
- `high15_bootstrap_top_n_frequency`: resample 중 High15 Top N에 포함된 비율
- `high15_technical_stability_pass`: 기본적으로 reference 우위 확률이 0.90 이상

이 bootstrap은 **NGS read sampling 안정성만** 반영합니다. PCR bias, cell sampling,
biological variation을 추정하지 않으며 biological confidence interval, p-value 또는
FDR로 표현하지 않습니다. PCR은 필요한 실험 단계로 인정하고, 낮은 cycle로 줄인
편향은 한계로 기록한 뒤 최종 후보를 개별 cloning으로 검증합니다.

## 후보 tier

| tier | 조건 | 사용법 |
|---|---|---|
| `tier1_clean_high_shift` | read support + High15가 reference 초과 + 기술적 안정성 + jackpot 아님 + weighted score도 reference 초과 | 우선 cloning |
| `tier2_high_tail` | 위 조건 중 weighted score support만 없음 | high-tail 후보로 cloning/분포 검토 |
| `tier3_review` | point estimate는 reference 초과하나 bootstrap 또는 jackpot 경고 존재 | 예비 후보, 낮은 우선순위 |
| `not_candidate` | reference 비초과 또는 eligibility 불충족 | 현재 cloning list 제외 |

`bin1`과 `bin2`가 각각 unsorted보다 농축되어야 한다는 조건은 후보 필터에 쓰지
않습니다. 매우 강한 UTR는 bin2를 지나 bin1에 집중되어 bin2가 오히려 감소할 수
있기 때문입니다. `top15_both_bins_enriched`는 진단 열로만 남습니다.

## 핵심 결과 파일

```text
results/sortseq/high15_primary_ranking.csv
results/sortseq/high15_candidates.csv
results/sortseq/top_candidates_for_cloning.csv
results/sortseq/top50_candidates_for_cloning.csv
results/sortseq/utr_results_full.tsv
```

- `high15_primary_ranking.csv`: eligibility 통과 UTR 전체의 robust High15 순위
- `high15_candidates.csv`: High15 point estimate가 comparator보다 높은 전체 후보
- `top_candidates_for_cloning.csv`: tier1·tier2만 모은 실제 cloning 우선순위
- `default_top_n_cloning_shortlist=TRUE`: 기본 robust-rank Top 50 중 tier1·tier2인 행
- `top50_candidates_for_cloning.csv`: 위 TRUE 행만 바로 연 파일(Top N 설정 변경 시 파일명도 변경)
- `top15_enrichment_ranking.csv`: 과거 workflow 호환용 별칭이며 내용은 새 High15 순위

R plotting 결과:

```text
results/sortseq/figures/top15/13_high15_bin1_bin2_structure.png
results/sortseq/figures/top15/14_high15_probability_vs_weighted_score.png
results/sortseq/figures/top15/15_top15_ranked_bin_probability_heatmap.png
results/sortseq/figures/top15/16_top15_reference_position.png
results/sortseq/figures/top15/top15_candidate_figures.pdf
```

## 기존 결과에서 다시 계산

Rescue와 LibraryQC는 다시 실행하지 않습니다.

```bash
# Python 환경
bash run_pipeline.sh reanalyze-analysis

# 별도 R 환경
bash run_pipeline.sh plot
```

기존 `results/sortseq`만 `archive/`로 이동되고 FASTQ, index rescue,
`results/library_qc`는 유지됩니다.

## 해석 한계와 후속 검증

Biological replicate가 없으므로 이 결과는 개별 construct 검증을 위한
prioritization입니다. `significant`, FDR, “translation rate가 증가했다”, 또는
fluorescence가 몇 배 증가했다는 결론에는 사용할 수 없습니다. 최종 판단은
상위 tier와 분포가 서로 다른 후보를 포함해 single-clone flow cytometry로 검증합니다.

## 분석 방식 참고

- Matreyek et al. (2018): 여러 fluorescence bin의 variant 분포를 weighted average로
  요약하는 방식. <https://pmc.ncbi.nlm.nih.gov/articles/PMC5980760/>
- Peterman and Levine (2016): threshold/high-tail enrichment와 mean fluorescence가
  서로 다른 phenotype 요약일 수 있음을 설명. <https://link.springer.com/article/10.1186/s12864-016-2533-5>
- Cao et al. (2021): pooled 5′UTR library에서 상위 fluorescence bin enrichment를
  후보 선별에 활용한 사례. <https://www.nature.com/articles/s41467-021-24436-7>
