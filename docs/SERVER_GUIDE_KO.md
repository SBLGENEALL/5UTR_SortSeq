# 서버 실행 가이드

## 1. 입력 구조

권장 배치는 다음과 같습니다.

```text
/data/user/MCET03/03_NGS/02_5UTR_sorting/
├── raw_data/                 # 원본: 절대 수정하지 않음
├── config/                   # 실제 index/경로 설정
├── 5UTR_SortSeq/             # 이 ZIP을 푼 폴더
├── results/                  # 자동 생성
└── archive/                  # --replace 때 이전 결과 보관
```

`raw_data`에는 이미 demultiplex된 bin1–6, unsorted R1/R2 FASTQ와 `Undetermined` R1/R2 FASTQ가 있어야 합니다. sample FASTQ가 7개라는 말은 보통 7개 library를 뜻하며, paired-end라면 실제 파일 수는 lane/chunk에 따라 14개 이상일 수 있습니다.

## 2. sample 이름 자동 추출

```bash
python scripts/init_project.py \
  --raw-dir /data/user/MCET03/03_NGS/02_5UTR_sorting/raw_data \
  --config-dir /data/user/MCET03/03_NGS/02_5UTR_sorting/config \
  --run-sample-sheet /data/user/MCET03/03_NGS/02_5UTR_sorting/raw_data/Analysis/1/Data/260812_sample_sheet.csv
```

다음 세 파일이 생깁니다.

- `detected_fastq_samples.csv`: 파일명에서 읽은 7개 sample prefix
- `SampleSheet.csv`: 원래 run SampleSheet에서 검증·추출한 실제 i7/i5
- `sample_map.csv`: NGS count matrix 열과 bin 의미 연결

`UTR_bin1_A2_s1_R1_001.fastq.gz`의 경우:

```csv
sample_id,ngs_column,sample_type,bin_number,population_fraction
bin1,UTR_bin1_A2,bin,1,0.05
```

처럼 연결됩니다. 파일명의 `A2`나 `s1` 자체가 dual index 서열이라는 뜻은 아닙니다. 실제 index는 `SampleSheet.csv`와 FASTQ header 마지막의 `I7+I5`를 사용합니다.

## 3. 실제 index 자동 추출

파이프라인은 원래 run SampleSheet의 `[Data]` 또는 `[BCLConvert_Data]` 구역에서
대소문자와 관계없이 `Sample_ID`, `index`, `index2` 열을 읽습니다. 각
`Sample_ID`는 FASTQ prefix와 정확히 또는 안전하게 정규화했을 때 1:1로
일치해야 합니다.

```csv
[Data]
Sample_ID,index,index2
UTR_bin1_A2,ACTUAL_I7,ACTUAL_I5
```

일치하지 않거나 중복이면 index를 추측하지 않고 setup이 중단됩니다.
`--run-sample-sheet`를 생략하면 `raw_data/Analysis` 아래의 `*sample_sheet*.csv`
또는 `*SampleSheet*.csv`가 정확히 하나일 때 자동 사용합니다. 후보가 여러 개면
명시 경로를 요구합니다.

`index2`의 방향은 주문서/샘플시트 표기와 FASTQ header 표기가 장비·conversion 설정에 따라 다를 수 있습니다. `I5_ORIENTATION=auto`가 as-given과 reverse-complement를 비교하고 선택합니다.

## 4. sorter 수치 입력

`sample_map.csv`의 nominal fraction은 현재 설계 `0.05, 0.10, 0.15, 0.20, 0.30, 0.20`으로 생성됩니다. sorter report에 실제 `collected events`가 있으면 `cells_collected` 열에 입력하세요. 여섯 bin 모두 값이 있을 때 파이프라인이 그 비율을 우선 사용합니다.

각 bin의 실제 median `log10(mCherry MFI)`가 있으면 `representative_log10_mfi`에 입력할 수 있습니다. 없으면 6→1 ordinal score를 사용합니다.

## 5. project env 확인

```bash
cp config/project.env.example \
  /data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env
```

특히 다음 두 값은 이전 plasmid run에서 실제 사용한 파일을 가리켜야 합니다.

```bash
LIBRARYQC_PY=/실제/NGS_LibraryQC/scripts/amplicon_qc_parallel.py
LIBRARYQC_CONFIG=/이전/plasmid/run/config/libraryqc.ini
```

reference, primer/anchor, exact/near-match 규칙을 plasmid와 cellular sample에서 동일하게 유지해야 enrichment bias를 줄일 수 있습니다.

## 6. preflight 판정

```bash
export SORTSEQ_PROJECT_CONFIG=/data/user/MCET03/03_NGS/02_5UTR_sorting/config/sortseq.env
bash run_pipeline.sh preflight
```

`results/index_rescue/rescue_inspection.json`에서 확인할 값:

- `missing_expected_samples_in_assigned_fastqs`: 비어 있어야 함
- `unexpected_assigned_fastq_samples`: 비어 있어야 함
- `headers_scanned_with_dual_index`: 0보다 커야 함
- `selected_i5_orientation`: 자동 선택 결과 기록
- `minimum_expected_dual_index_distance`: rescue radius와 비교

header에 `i7+i5`가 없으면 FASTQ만으로 read-level rescue할 수 없습니다. 이 경우 synchronized I1/I2 FASTQ 또는 original BCL에서 다시 demultiplex해야 합니다.

## 7. 계산 리소스 확인

```bash
bash run_pipeline.sh resources
```

이 서버에서는 `physical_cores: 128`, `logical_cpus: 256`, `numa_nodes: 4`,
`libraryqc_workers: 128`이 기대값입니다. `rescue_compression`은 `pigz`가 이미 설치된
서버에서는 `pigz`, 그렇지 않으면 `parallel_python`으로 표시됩니다. 두 경우 모두
오프라인에서 추가 package 설치 없이 실행할 수 있습니다.

## 8. 실행과 재실행

각 단계의 결과를 이해하고 확인하기 쉽도록 다음처럼 하나씩 실행하는 것을
권장합니다. 자세한 목적과 완료 판정은
[PIPELINE_WORKFLOW_KO.md](PIPELINE_WORKFLOW_KO.md)를 참고하세요.

```bash
bash run_pipeline.sh preflight
bash run_pipeline.sh rescue
bash run_pipeline.sh libraryqc
bash run_pipeline.sh analyze
bash run_pipeline.sh plot
```

완료된 단계는 다시 실행하지 않고 다음 단계부터 이어갑니다. 기존 결과가 있으면
기본적으로 중단합니다. 모든 결과를 archive로 보존하고 처음부터 다시 실행할 때만:

```bash
bash run_pipeline.sh full --replace
```

## 9. 실패 시 빠른 확인

| 오류 | 먼저 볼 것 |
|---|---|
| SampleSheet sample mismatch | `detected_fastq_samples.csv`와 `Sample_ID`가 같은지 |
| No i7+i5 strings | FASTQ header 두 번째 token 마지막 필드 |
| matrix column not found | `sample_map.csv`의 `ngs_column`과 matrix header |
| R package missing | 인터넷 가능한 별도 환경에서 준비하거나 기존 pDNA R 환경 확인 |
| coverage passing UTR 없음 | NGS_LibraryQC assigned reads와 cutoff 50/100 |
