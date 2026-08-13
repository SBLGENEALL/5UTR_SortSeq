# Undetermined dual-index rescue 원칙

큰 `Undetermined` FASTQ를 통째로 버리지는 않지만, 모든 read를 억지로 7개 sample에 넣지도 않습니다.

## 배정 규칙

관측된 FASTQ header index `i7+i5`와 7개 expected pair의 Hamming distance를 계산합니다. 다음 조건을 모두 만족할 때만 rescue합니다.

1. 가장 가까운 expected pair가 하나뿐이다.
2. i7+i5 총 mismatch가 기본 2 이하이다.
3. 각 index mismatch가 기본 2 이하이다.
4. 두 번째 후보와의 거리 차이가 `MIN_INDEX_DISTANCE_MARGIN` 이상이다.
5. paired-end이면 R1/R2 read ID와 header index가 서로 일치한다.

나머지는 `Undetermined_residual`로 남습니다. residual은 NGS_LibraryQC QC에는 나타날 수 있지만 bin1–6/unsorted phenotype 분석에는 사용하지 않습니다.

## 원본 보존

- 기존 7개 assigned FASTQ는 결과 폴더에 symlink합니다.
- rescue read만 별도 `L900` chunk로 새로 씁니다.
- 원본 `raw_data`는 수정하거나 덮어쓰지 않습니다.
- 동률이나 index 길이 불일치는 residual로 보존합니다.

## 2 mismatch 사용 시 QC

Expected dual-index pair 사이의 최소 거리가 rescue radius의 두 배 이하이면 서로의 error sphere가 가까워집니다. 파이프라인은 경고를 남기며, unique-nearest 규칙은 그대로 유지합니다.

최종 판단에서는 다음을 확인하세요.

- rescue가 한 sample에 비정상적으로 몰리지 않는가
- `rescued_distance_2`가 0/1 mismatch보다 과도하게 많지 않은가
- `ambiguous_nearest_index`가 많은가
- 상위 unknown barcode가 실제 expected index의 1–2 mismatch로 설명되는가
- 2 mismatch를 포함했을 때와 1 mismatch까지만 포함했을 때 hit tier가 크게 바뀌지 않는가

민감도 분석이 필요하면 private env에서 `MAX_TOTAL_INDEX_MISMATCHES=1`로 별도 결과 디렉터리에 한 번 더 실행해 비교합니다.
