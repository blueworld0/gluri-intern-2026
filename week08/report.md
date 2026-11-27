# Week08 다른 날짜·궤도 검증

4~5월 후보 조사와 고도·경사 비교는 [week08 후속 분석](spring_report.md)에 정리했다.

2025년 4월 25일 GEDI 산림 184점으로 week07의 10밴드 랜덤 포레스트를 다시 적합하고, 같은 연구 구역의 2025년 2월 4일 GEDI(궤도 O34820) 산림 35점에 적용했다. 같은 35점에서 랜덤 포레스트의 RMSE는 **46.35 Mg/ha**, 4월 학습 평균 예측의 RMSE는 **42.77 Mg/ha**였다. 고 AGBD(100 Mg/ha 이상) 10점은 두 모델 모두 과소 예측했다. [예측 결과](outputs/february_holdout_predictions.csv) · [상세 지표](outputs/february_holdout_summary.json)

| 평가 범위 | 지표 | 랜덤 포레스트 | 학습 평균 |
| --- | --- | ---: | ---: |
| 전체 35점 | RMSE (Mg/ha) | 46.35 | 42.77 |
| 전체 35점 | MAE (Mg/ha) | 38.15 | 33.83 |
| 전체 35점 | R² | −0.179 | −0.004 |
| 고 AGBD 10점 | 평균 과소 예측 (Mg/ha) | 62.34 | 57.38 |

## 표본과 영상

- GEDI 품질 통과 후보 533점 중 week07과 동일한 EGIS 중심점 분류로 산림 220점을 확인했다. 이 중 61점이 고 AGBD였다. [산림 분류](outputs/february_forest_summary.json)
- 1월 28일~2월 11일 UTC, 경도 128.75~129.00°E 및 위도 35.00~35.25°N의 Earth Search 제품 12개에서 각 후보의 20m SCL을 확인했다. 2월 2일과 7일에는 타일 내 후보가 모두 구름·그림자 코드였다. 1월 28일 T52SDE는 220점을 모두 덮었고, SCL 식생 코드 4를 통과한 표본이 단일 제품 중 가장 많은 60점이었다. 고 AGBD 수나 모델 성능은 영상 선택에 쓰지 않았다. GEDI와 영상의 시차는 약 167시간이었다. [제품별 집계](outputs/february_scl_products.csv) · [Earth Search 원응답](outputs/february_earth_search.json)
- 원본 JP2의 SCL은 60점의 중심에서 공개 SCL GeoTIFF와 일치했다. 원형 영역(반경 12.5m) 전체에서 SCL 4·5·6, CLASSI 3개 층이 모두 0, 10개 밴드 유효 반사율 조건을 적용했다. 60점 중 25점은 원의 일부가 CLASSI를 통과하지 않아 제외됐다. 최종 35점 중 10점이 고 AGBD였다. [점별 품질 기록](outputs/february_footprint_audit.csv)
- 4월과 2월의 평가·학습 ID는 겹치지 않았고, 최종 평가점에서 가장 가까운 학습점까지 최소 8.26km였다. 2월 GEDI와 1월 말 영상으로 구성된 단일 평가 자료이므로 관측일·궤도·계절·위치 영향을 따로 추정할 수 없다. 비교 기준은 GEDI 추정 AGBD다.

## 원자료와 재현

공식 제품명은 `S2C_MSIL2A_20250128T021001_N0511_R103_T52SDE_20250128T060411.SAFE`이며 CDSE UUID는 `7b7f1032-3e93-479a-9d73-bb922ce81888`이다. [CDSE 제품 메타데이터](outputs/february_selected_cdse.json)

SCL GeoTIFF는 Earth Search 원응답의 각 `assets.scl.href`에서 내려받아 `source/scl/<item_id>.tif`에 두었다. 원본 JP2는 공개 AWS 버킷 `sentinel-s2-l2a`의 `tiles/52/S/DE/2025/1/28/0/` 아래에서 받았다. `R10m`의 B02·B03·B04·B08, `R20m`의 B05·B06·B07·B8A·B11·B12·SCL, `qi`의 CLASSI_B00 및 밴드별 QUALIT·DETFOO를 `source/january_selected/`에 파일명으로 저장했다. 제품 경로는 타일의 `productInfo.json`에서 확인했고, 해당 경로의 `metadata.xml`을 `product_metadata.xml`로 저장했다. 원본 체크섬은 상세 지표 JSON에 기록했다. [공개 AWS 자료 안내](https://registry.opendata.aws/sentinel-2/) · [Earth Search 자료 안내](https://registry.opendata.aws/sentinel-2-l2a-cogs/)

저장소 루트에서 실행한다. 재현 과정에서 원본과 기존 결과를 보존하려면 두 실행 모두 별도 출력 폴더를 지정한다.

```powershell
& .\week04\.venv\Scripts\python.exe -B .\week08\screen_february_sentinel.py --output-dir .\week08\scl_reproduction
& .\week04\.venv\Scripts\python.exe -B .\week08\evaluate_february_holdout.py --scl-points .\week08\scl_reproduction\february_scl_points.csv --output-dir .\week08\holdout_reproduction
```
