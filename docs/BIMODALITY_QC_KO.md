# UTR A형 양극화 분포 QC

## 먼저 결론

bin1·2와 bin5·6에 동시에 많고 bin3·4가 비어 보이는 UTR를 곧바로 “좋은 UTR” 또는
“두 개의 biological state를 만든 UTR”라고 부르지 않습니다. 현재 6-bin pooled
Sort-seq에서 말할 수 있는 것은 **read-supported high+low-tail 분포 후보**입니다.

- 목표가 일부 cell의 매우 높은 mCherry라면: 후속 검증할 가치가 있는 high-tail 후보
- 목표가 모든 cell에서 균일하고 높은 mCherry라면: 낮은 cell도 많으므로 우선 후보가 아님
- 총 read가 낮을 때만 보이면: sampling/PCR jackpot 가능성이 큼
- 높은 read에서도 반복되고 개별 flow cytometry에서 두 peak가 보이면: biological
  heterogeneity 후보로 해석 가능

6개 구간의 pooled read만으로 연속 fluorescence 분포의 peak 두 개를 직접 증명할 수
없으므로 결과에서는 `bimodality-like`, `high+low-tail`, `polarized shape`라고 표현합니다.

## 계산하는 값

기존 분석에서 depth와 sorter-bin size를 보정한 UTR별 bin probability를 사용합니다.

```text
high tail   H = p_bin1 + p_bin2
middle      M = p_bin3 + p_bin4
low tail    L = p_bin5 + p_bin6

extreme polarization index = 2 × min(H, L)
```

Polarization index는 양쪽 tail 중 더 작은 쪽까지 얼마나 충분히 존재하는지 나타냅니다.
한쪽에만 몰리면 작고, high와 low가 각각 50%이면 최대 1입니다. 이 지수만으로 후보를
정하지 않고 middle valley와 raw-read support를 함께 봅니다.

기본 `clear_bimodal_flag=TRUE` 기준은 다음과 같습니다.

```text
total six-bin raw count >= 200
bin1+bin2 raw count >= 20
bin5+bin6 raw count >= 20
H >= 0.20
L >= 0.20
M <= 0.30
max(p3,p4) / min(max(p1,p2), max(p5,p6)) <= 0.75
```

`utra_like_strong_polarization_flag`는 더 UTR A와 비슷한 모양을 찾는 보조 flag입니다.

```text
H >= 0.20, L >= 0.50, M <= 0.20
위 raw-read support와 valley 조건 통과
```

이 기준은 통계적 유의성이나 FDR이 아니라 현재 데이터의 분포 모양과 read-count
의존성을 검토하기 위한 exploratory QC cutoff입니다.

## 실행 순서

기존 rescue, LibraryQC, analyze 결과를 다시 만들지 않습니다. Python 환경에서 다음
추가 QC만 실행합니다.

```bash
cd /data/user/MCET03/03_NGS/02_5UTR_sorting/5UTR_SortSeq-main
export SORTSEQ_PROJECT_CONFIG=/data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env

bash run_pipeline.sh bimodality-qc 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/bimodality_qc.log
```

그다음 별도 R 환경에서 그림만 다시 생성합니다.

```bash
bash run_pipeline.sh plot 2>&1 | tee \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/plot_bimodality.log
```

## 결과 파일

Python QC:

```text
results/sortseq/bimodality_qc/bimodality_summary.csv
results/sortseq/bimodality_qc/bimodality_all_variants.csv
results/sortseq/bimodality_qc/clear_bimodal_candidates.csv
results/sortseq/bimodality_qc/utra_like_strong_polarization.csv
results/sortseq/bimodality_qc/bimodality_by_read_count.csv
results/sortseq/bimodality_qc/bimodality_count_sensitivity.csv
```

R 그림:

```text
results/sortseq/figures/bimodality/17_polarization_vs_read_count.png
results/sortseq/figures/bimodality/18_bimodal_fraction_by_read_count.png
results/sortseq/figures/bimodality/19_high_vs_low_tail_probability.png
results/sortseq/figures/bimodality/20_clear_bimodal_bin_probability_heatmap.png
results/sortseq/figures/bimodality/bimodality_qc_figures.pdf
```

## read count 경향을 판정하는 법

다음 두 표를 같이 봅니다.

```bash
column -s, -t results/sortseq/bimodality_qc/bimodality_by_read_count.csv
column -s, -t results/sortseq/bimodality_qc/bimodality_count_sensitivity.csv
```

| 결과 | 해석 |
|---|---|
| `200–499`에서 많고 `>=1000`에서 급감 | low-count sampling/jackpot 가능성이 큼 |
| 최소 count를 200→500→1000으로 높여도 비율이 비슷 | read count만으로 설명하기 어려움 |
| count와 polarization의 Spearman rho가 뚜렷한 음수 | 낮은 count일수록 양극화가 커지는 경향 |
| 높은 count에서도 clear 후보가 남고 replicate에서 재현 | biological heterogeneity 가설이 강해짐 |

Spearman p-value가 작아도 biological replicate가 하나이면 독립적인 생물학적 재현성을
뜻하지 않습니다. UTR A형 후보는 개별 construct로 다시 transfection한 뒤 single-cell
flow cytometry histogram, replicate, 필요하면 mRNA abundance를 확인합니다.

## 문헌과의 연결

- Matreyek et al. (2018)은 여러 fluorescence bin의 variant abundance를 weighted
  average로 요약했습니다. <https://pmc.ncbi.nlm.nih.gov/articles/PMC5980760/>
- Peterman and Levine (2016)은 threshold enrichment와 mean fluorescence가
  cell-to-cell variability에 따라 서로 다른 후보를 고를 수 있음을 설명했습니다.
  <https://link.springer.com/article/10.1186/s12864-016-2533-5>
- Gilliot and Gorochowski (2023)는 Flow-seq의 낮은 read/cell coverage에서 phenotype
  추정과 mode-bin 분류가 불안정해질 수 있음을 평가했습니다.
  <https://academic.oup.com/bioinformatics/article/39/5/btad277/7135834>

