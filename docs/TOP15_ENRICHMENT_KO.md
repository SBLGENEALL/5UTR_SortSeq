# bin1·bin2 중심 고발현 5′UTR 랭킹

## 왜 별도 랭킹을 사용하는가

`expected_bin_score`는 여섯 bin의 전체 분포를 한 숫자로 요약합니다. 따라서 bin3에
많이 존재하는 안정적인 intermediate-high UTR도 높은 점수를 받을 수 있습니다.
이번 실험의 hit-selection 목표는 평균적인 위치보다 **가장 높은 mCherry 구간인
bin1+bin2에 unsorted보다 농축된 UTR**을 찾는 것입니다.

기존 `expected_bin_score`를 삭제하거나 bin1·bin2에 임의의 큰 숫자를 주지 않습니다.
두 지표를 다음처럼 역할별로 분리합니다.

| 지표 | 답하는 질문 |
|---|---|
| `expected_bin_score` | target gate 안에서 전체적으로 어느 fluorescence 위치에 있는가? |
| `top15_vs_unsorted_log2_enrichment` | whole unsorted와 비교해 bin1+bin2에 얼마나 농축됐는가? |

## 계산식

UTR `i`의 sample별 raw count를 `c`, 해당 sample 전체 assigned UTR read를 `N`이라
하면 depth-normalized frequency는 다음과 같습니다.

```text
f_i1 = c_i1 / N_1
f_i2 = c_i2 / N_2
u_i  = c_iu / N_u
```

bin1과 bin2의 sorter fraction을 각각 `w1`, `w2`라 하면 combined top-bin frequency는:

```text
f_high15 = (w1 × f_i1 + w2 × f_i2) / (w1 + w2)
```

현재 nominal fraction에서는 `w1=0.05`, `w2=0.10`, `w1+w2=0.15`입니다. 최종
unsorted 대비 enrichment와 log2 enrichment는:

```text
top15_vs_unsorted_enrichment      = f_high15 / u_i
top15_vs_unsorted_log2_enrichment = log2(f_high15 / u_i)
```

0 count에서 무한값이 생기지 않도록 sample frequency에는 기본 pseudocount 0.5를
적용합니다.

| log2 enrichment | 해석 |
|---:|---|
| 0 | unsorted와 같은 representation |
| 0.585 | 1.5배 농축 |
| 1 | 2배 농축 |
| 2 | 4배 농축 |

## 후보와 priority 후보

기본 read-support 조건은 다음과 같습니다.

```text
unsorted count >= 50
total six-bin count >= 200
bin1 + bin2 raw count >= 20
```

`top15_candidate_flag=TRUE`가 되려면 read support를 통과하고, bin1과 bin2가 각각
unsorted 대비 1배 이상이며, combined top15 enrichment가 original보다 커야 합니다.
Reference가 없으면 pool-neutral 값인 1을 comparator로 사용합니다.

`top15_priority_candidate_flag=TRUE`는 위 조건에 더해
`single_bin_jackpot_suspect=FALSE`인 후보입니다. Jackpot flag가 있는 UTR도 결과에서
삭제하지 않고 review 대상으로 남깁니다.

Biological replicate가 하나뿐이므로 이 flag는 통계적으로 유의한 hit 또는 FDR을
뜻하지 않습니다. 개별 construct flow cytometry 검증을 위한 exploratory rank입니다.

## 결과 파일

Python 분석 후 다음 UTF-8 CSV가 생성됩니다.

```text
results/sortseq/top15_enrichment_ranking.csv
results/sortseq/top15_candidates.csv
results/sortseq/top15_priority_candidates.csv
```

- `top15_enrichment_ranking.csv`: read-support를 통과한 전체 UTR를 top15 enrichment 순 정렬
- `top15_candidates.csv`: bin1·2 일관성 및 comparator 조건을 통과한 후보
- `top15_priority_candidates.csv`: jackpot suspect를 제외한 우선 검증 후보

R plotting 후에는 다음 파일이 추가됩니다.

```text
results/sortseq/figures/top15/13_bin1_bin2_unsorted_enrichment.png
results/sortseq/figures/top15/14_top15_enrichment_vs_expected_score.png
results/sortseq/figures/top15/15_top15_ranked_bin_probability_heatmap.png
results/sortseq/figures/top15/16_top15_reference_position.png
results/sortseq/figures/top15/top15_candidate_figures.pdf
```

## 기존 결과에서 다시 계산

Rescue와 LibraryQC는 다시 실행하지 않습니다. Python 환경과 R 환경이 분리되어 있으므로
다음처럼 나누어 실행합니다.

```bash
# Python 환경: 기존 sortseq 결과만 archive하고 step 4 재계산
bash run_pipeline.sh reanalyze-analysis

# R 환경: 새 표에서 step 5 그림 생성
bash run_pipeline.sh plot
```

## 참고한 분석 방식

- Cao et al. (2021)은 5′UTR library를 상위 GFP bin들로 sorting한 뒤 각 top bin을
  unsorted에 대한 log2 enrichment로 비교하고, 모든 top bin에서 일관되게 농축된
  후보를 검증했습니다.
  <https://www.nature.com/articles/s41467-021-24436-7>
- Matreyek et al. (2018)의 VAMP-seq 방식은 모든 bin 분포의 weighted average를 사용해
  평균 phenotype을 추정합니다. 이는 현재의 `expected_bin_score`에 해당하는 보조
  관점입니다. <https://pmc.ncbi.nlm.nih.gov/articles/PMC5980760/>
- Peterman and Levine (2016)은 top threshold enrichment와 mean fluorescence가 서로
  다른 요약이며, cell-to-cell variability에 따라 관계가 달라질 수 있음을 설명합니다.
  <https://link.springer.com/article/10.1186/s12864-016-2533-5>

