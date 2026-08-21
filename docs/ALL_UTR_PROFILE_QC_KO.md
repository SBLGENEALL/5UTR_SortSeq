# 전체 UTR normalized-enrichment profile 시각화

## 무엇을 보여주는가

이 QC는 모든 UTR에 대해 다음 값을 계산해 CSV로 저장합니다.

```text
q_ib = f_ib / sum_k(f_ik)
     = (p_ib/w_b) / sum_k(p_ik/w_k)
```

`q`는 UTR마다 bin1–6 합이 1인 **normalized relative-enrichment shape**입니다.
`p=P(bin|UTR, gate)`와 달리 실제 세포 probability로 부르지 않습니다.

전체 profile clustering과 overview 그림에 사용하는 기준은 다음 하나입니다.

```text
6-bin raw count 합 >= 200
```

정확히 200인 UTR도 포함됩니다. 이 전체-profile overview에는 `bin1+bin2`나 unsorted
cutoff를 적용하지 않습니다. 그래야 high candidate뿐 아니라 bin4–6 쪽으로 이동한
low-expression profile도 함께 볼 수 있습니다. 반면 Step 6·7·8 후보 순위 비교는
high-tail 안정성을 위해 별도로 `bin1+bin2 >= 20`을 적용합니다.

각 UTR의 여섯 `q` 합은 1입니다. 예를 들어 다음 세 UTR은 서로 다른 normalized
enrichment 모양을 가집니다.

| UTR | bin1 | bin2 | bin3 | bin4 | bin5 | bin6 |
|---|---:|---:|---:|---:|---:|---:|
| A | 20% | 20% | 20% | 20% | 20% | 0% |
| B | 50% | 10% | 10% | 10% | 10% | 10% |
| C | 40% | 30% | 20% | 10% | 0% | 0% |

## 실행 순서

Python 환경과 R 환경을 나누어 실행합니다.

```bash
export SORTSEQ_PROJECT_CONFIG=/data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env

# Python 환경
bash run_pipeline.sh profile-qc 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/all_utr_profile_qc.log

# R 환경으로 전환한 뒤
bash run_pipeline.sh plot 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/all_utr_profile_plot.log
```

`analyze`를 새로 실행하면 `profile-qc`도 자동 실행됩니다. 기존 v0.2.0 결과가 이미
있다면 위 두 명령만 실행하면 되며 rescue, LibraryQC, 기존 scoring 결과는 지워지지
않습니다.

## Python 결과표

```text
results/sortseq/all_utr_profiles/all_utr_profile_assignments.csv
results/sortseq/all_utr_profiles/all_utr_profile_cluster_summary.csv
results/sortseq/all_utr_profiles/all_utr_normalized_enrichment_profiles.csv
results/sortseq/all_utr_profiles/top200_normalized_enrichment_profiles.csv
results/sortseq/all_utr_profiles/top_normalized_enrichment_profiles.csv
results/sortseq/all_utr_profiles/all_utr_profile_manifest.json
```

- `all_utr_normalized_enrichment_profiles.csv`: 약 2,001개 **모든 UTR**의 q 값과
  coverage/profile eligibility, High15 rank, original 표시
- `top200_normalized_enrichment_profiles.csv`: High15 상위 200개와 `original`
  (`original`이 Top 200 밖이어도 비교용으로 마지막에 추가)

### `all_utr_profile_assignments.csv`

한 행이 한 UTR입니다.

| 열 | 의미 |
|---|---|
| `profile_cluster` | 비슷한 6-bin 분포끼리 묶은 cluster 번호 |
| `profile_label` | cluster 번호와 포함 UTR 수 |
| `heatmap_order` | heatmap에서의 행 순서 |
| `bin1_probability`–`bin6_probability` | UTR별 추정 bin 분포; 합계 1 |
| `bin1_normalized_enrichment_share`–`bin6_normalized_enrichment_share` | `q=f/sum(f)`; 합계 1 |
| `equal_bin_high_share` | q1+q2; Top15 probability가 아닌 profile 보조값 |
| `high15_probability` | bin1 + bin2 probability |
| `expected_bin_score` | 6–1 weighted score |
| `is_reference_variant` | original/orginal 여부 |

이 CSV를 열면 각 cluster에 어떤 UTR이 들어갔는지 정확히 확인할 수 있습니다.

### `all_utr_profile_cluster_summary.csv`

한 행이 한 profile cluster입니다. `binN_mean_normalized_enrichment_share`는 그 cluster
UTR들의 평균 q profile이고 `utr_count`는 포함 UTR 수입니다.
`reference_present=TRUE`인 행이 original이 속한 cluster입니다.

## R 그림

```text
results/sortseq/figures/all_utr_profiles/25_all_utr_bin_probability_heatmap.png
results/sortseq/figures/all_utr_profiles/26_all_utr_relative_enrichment_heatmap.png
results/sortseq/figures/all_utr_profiles/27_profile_cluster_composition.png
results/sortseq/figures/all_utr_profiles/28_profile_cluster_variability.png
results/sortseq/figures/all_utr_profiles/29_all_utr_normalized_enrichment_profile_heatmap.png
results/sortseq/figures/all_utr_profiles/30_normalized_enrichment_profile_cluster_means.png
results/sortseq/figures/all_utr_profiles/31_normalized_enrichment_profile_cluster_variability.png
results/sortseq/figures/all_utr_profiles/32_top200_normalized_enrichment_profile_heatmap.png
results/sortseq/figures/all_utr_profiles/33_top200_normalized_enrichment_profiles_page_01.png
results/sortseq/figures/all_utr_profiles/top200_normalized_enrichment_profile_pages.pdf
results/sortseq/figures/all_utr_profiles/all_utr_profile_figures.pdf
```

### 25번: 절대 probability heatmap

- 가로축: bin1–bin6
- 세로축: UTR; 한 줄이 UTR 하나
- 색: 해당 UTR의 추정 bin probability
- 왼쪽 구획: 비슷한 분포끼리 묶인 profile cluster
- 빨간 테두리: original/orginal

UTR 이름 1,000개 이상을 세로축에 모두 쓰면 읽을 수 없으므로 이름은 숨깁니다. 정확한
UTR 이름과 행 순서는 assignments CSV의 `heatmap_order`에서 찾습니다.

### 26번: bin 기본 크기 대비 상대 농축 heatmap

표시값은 다음과 같습니다.

```text
log2[ P(bin b | UTR, target gate) / sorter bin fraction_b ]
```

- 빨강: 해당 bin의 기본 크기보다 UTR가 더 많이 존재
- 흰색: bin 기본 크기와 비슷
- 파랑: bin 기본 크기보다 적게 존재

예를 들어 bin5 probability가 20%이면 절대 heatmap에서는 커 보일 수 있지만,
bin5의 기본 크기가 30%이므로 상대 농축값은 `log2(0.20/0.30)=-0.585`이고 파란색입니다.

### 27번: cluster 평균 구성

각 profile cluster의 평균 6-bin 분포를 100% 누적 막대로 요약합니다. 어떤 분포 유형이
전체 library에서 큰 비중을 차지하는지 빠르게 확인할 수 있습니다.

### 28번: cluster 내부 변이

- 가는 회색선: UTR 하나의 여섯-bin 분포
- 굵은 주황선: cluster 평균
- 빨간선: original/orginal

평균만 보면 숨겨질 수 있는 cluster 내부의 넓은 변이와 outlier를 확인합니다.

### 29–31번: normalized enrichment profile

- 29번: `q=f/sum(f)` 전체 UTR heatmap
- 30번: 8개 profile cluster의 평균 q 구성
- 31번: cluster 안의 개별 q profile과 평균선
- 빨간 테두리·별·선: `original` 위치

### 32–33번: High15 Top 200

- 32번: High15 상위 200개의 q heatmap; `bin1`이 왼쪽
- `original`이 Top 200 밖이어도 마지막 행에 추가하고 빨간색으로 표시
- 33번: 25개 후보씩 나눈 개별 profile PNG
- `top200_normalized_enrichment_profile_pages.pdf`: 모든 page를 합친 PDF
- 각 후보 panel의 빨간 점선은 `original`, 주황선은 후보입니다.

## 군집 방법과 주의점

여섯 q 값에 제곱근 변환을 적용한 뒤 Hellinger distance에 해당하는 공간에서
Ward hierarchical clustering을 수행합니다. 기본 cluster 수는 8이며 설정에서 바꿀 수
있습니다.

```bash
ALL_UTR_PROFILE_CLUSTERS=8
ALL_UTR_PROFILE_MIN_TOTAL_COUNT=200
ALL_UTR_PROFILE_TOP_N=200
```

Profile 1은 cluster median `q1+q2`가 가장 높은 쪽이 되도록 번호만 재배열합니다.
따라서 **profile cluster 번호는 후보 rank나 통계적으로 발견된 biological class가
아닙니다.** 전체 분포를 탐색하기 위한 요약입니다. 최종 cloning 순위는 계속
`high15_final_rank`와 `top_candidates_for_cloning.csv`를 사용합니다.
