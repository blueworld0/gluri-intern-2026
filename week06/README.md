# week06 · 환산 계수를 국가 고유계수로 확정 (금정산 마무리)

5주차까지 임시값이던 **목재기본밀도(WD)·바이오매스 확장계수(BEF)·뿌리함량비(R)** 를
출처가 있는 **국가 고유계수**(이선정 외 2015, 국립산림과학원 연구진, 표 1)로 바꾸고,
금정산 임상도의 수종 면적으로 가중평균해서 모든 숫자를 다시 계산했다.

## 실행 방법
- `week06/` 폴더를 `week04/`, `week05/` 폴더 **옆에** 둔다 (위성·임상도 자료는 week04/data 를 그대로 씀)
```bash
cd week06/src
python 06_national_coef.py     # 수종별 계수 + 가중평균 + 새 R 지도 (임상도 읽느라 1~2분)
python 07_carbon_national.py   # 위성 D 저장량·흡수량 다시 계산 (수십 초)
python 08_uncertainty.py       # 몬테카를로 불확도 2만 회 (약 10초)
python 09_stats_compare.py     # 국가 통계 vs 위성 (수 초)
python 10_loss_check.py        # 2022년 감소가 실제 숲 손실인지 확인 (수 초)
```

## 계수 출처
| 항목 | 값 |
|---|---|
| 논문 | 이선정·임종수·손영모·김래현 (2015). 강원도 산림의 임상별, 수종별 탄소저장량 및 온실가스 흡수량 산정. 한국기후변화학회지 6(4): 303–310. https://doi.org/10.15531/KSCCR.2015.6.4.303 |
| 표 | Table 1 "National specific expansion factors" (수종 11개의 WD, BEF, R) |
| 성격 | 저자 모두 국립산림과학원. 본문: 온실가스종합정보센터(GIR)의 검증을 받아 공표된 국가 고유계수 |

## 수종 → 계수 연결 규칙 (`src/config.py`)
- 표에 있는 수종(곰솔, 소나무, 리기다, 잣나무, 낙엽송, 상수리, 신갈, 굴참)은 그 값
- 기타참나무류 = 표의 참나무 4종(굴참·상수리·신갈·졸참) 평균
- 침활혼효림 = 곰솔 1/2 + 참나무 4종 1/2 (논문도 혼효림은 침엽·활엽 반반으로 계산)
- 편백·삼나무 = '기타 침엽수', 아까시·서어 등 = '기타 활엽수'

## 결과 요약 (임시 경계, 숲 6,623 ha)
| 항목 | 5주차 | **6주차 (확정)** |
|---|---|---|
| 평균 R | 0.321 | **0.327** |
| 1 m³ → 지상부 바이오매스 | 0.847 t (임시) | **0.874 t** (낮게~높게 0.789–0.897) |
| 1 m³ → 탄소 | 0.561 tC (임시, 범위 0.37–0.75) | **0.583 tC** (범위 0.516–0.606) |
| 저장량 2022 | 60.0만 tC | **60.2만 tC** (90% 42.4만–75.2만) |
| 연간 흡수량 | 1.37만 tC/년 | **1.38만 tC/년 = 5.0만 tCO₂/년** (90% 3.0만–7.1만 tCO₂) |
| 통계(전국) 흡수 속도 2015→2020 | 2.16 tC/ha/년 | **2.24 tC/ha/년** (1.98–2.33) vs 위성 2.08 |
| 2020 바이오매스 수준 (위성 139.7 t/ha) | 전국 139.9 | 전국 144.4 / **부산 156.9 / 금정구 154.7 / 북구 164.0** → 평균 방식으론 위성이 10–15% 낮지만, **임상별 보정 시 부산 −4.9%, 금정구 −3.9%** |
| 금정구 통계 흡수 속도 2015→2020 | — | 평균 방식 2.89 / 임상 보정 1.55 tC/ha/년 → 위성 2.08 이 사이에 있음 |
| 2022년 숲 손실 (Hansen) | — | 0.7 ha (평년 1.1 ha) → 2022 감소는 숲이 사라져서가 아님 |

## 결과 파일 (`outputs/`)
| 파일 | 내용 |
|---|---|
| `national_coef_source_table.csv` | 논문 표 1 값 그대로 (11개 수종) |
| `species_coef_table.csv` | 금정산 수종별 면적과 붙인 계수 |
| `conversion_factors_national.csv` | 부피 → 바이오매스·탄소 환산계수 (기본/낮게/높게/5주차 임시값) |
| `national_coef_summary.json` | 가중평균 WD·BEF·R 요약 |
| `stock_by_year_national.csv` | 연도별 저장량 (week04 대비 차이 포함) |
| `sink_national.csv` | 연간 흡수량 (week04 대비 차이 포함) |
| `uncertainty_summary.csv`, `uncertainty_contribution.csv`, `map_error_by_block.csv`, `uncertainty_settings.json`, `sink_uncertainty_hist.png` | 불확도 |
| `stats_vs_satellite_rate.csv`, `stats_vs_satellite_level.csv`, `stats_vs_satellite.png` | 국가 통계 비교 |
| `stats_vs_satellite_by_forest_type.csv`, `forest_type_factors.json` | 임상별(침엽·활엽·혼효) 보정 비교 |
| `stats_vs_satellite_rate_by_forest_type.csv` | 부산·구별 흡수 속도 임상 보정 (2015→2020) |
| `forest_loss_by_year.csv` | 경계 안 연도별 숲 손실 면적 (Hansen) |

## 자료
- `data_stats/forest_stock.csv`: ha당 임목축적. 전국 2010·2015·2020 + 부산광역시·금정구 2015·2020 + 북구 2020 (산림청 산림기본통계, KOSIS 「행정구역별 면적 및 축적」)
- `data_stats/busan_forest_type_2015_kosis.csv`, `busan_forest_type_2020_kosis.csv`: KOSIS 「행정구역별 임상별 면적 및 축적」 2015·2020 (부산, 내려받은 그대로)
- `data_extra/lossyear.tif`: Hansen Global Forest Change v1.12 손실 연도 (30 m, 금정산 주변만 잘라 둠, CC BY 4.0). 경도 129.00° 부터라 경계 서쪽 끝 약 4.5%는 빠짐
- 새 R 지도는 `week04/data/imsang/r_map_2019_national.tif` 로 저장됨 (깃허브에는 안 올림)
