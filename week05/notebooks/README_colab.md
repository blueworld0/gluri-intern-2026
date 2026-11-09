# Colab 실행 안내

파일: `NEON_DeepForest_AGB_2021.ipynb` (46셀: 코드 22, 설명 24). 나무 찾기 개선 A/B/C를 추가한 수정본입니다.

## 추가된 나무 찾기 비교

- 동일 23개 조사구 ID를 설정에 고정합니다. GRSM 15개 + ORNL 8개이며 누락·대체가 생기면 중단합니다.
- A: score 0.1/0.2/0.3/0.4/0.5 × patch 300/400/600. 바깥 조사구 GroupKFold 평가 정답은 조합 선택에 쓰지 않습니다.
  안쪽 GroupKFold로 보이는 나무 F1 최고 조합을 선택하고 바깥 OOF 탐지만 합쳐 성능표에 씁니다.
- 전 조사구에서 선택한 단일 조합은 새 자료 적용용입니다. 그 선택 점수를 독립 성능으로 보고하지 않습니다.
- B: 1m CHM Gaussian sigma=1m, 원형 가변 창 지름=2.51+0.00901H², 3m 이상 최대점, 시드 watershed.
  H는 탐지한 꼭대기 픽셀의 원본 CHM 높이, CD는 watershed 수관 면적의 등가원 지름입니다.
- C: A OOF 상자 안 최고 꼭대기로 H를 보완하고 상자 CD는 유지합니다. 꼭대기 없는 상자는 보존합니다.
  상자 밖 꼭대기는 B 수관으로 추가하고 1m 중복을 제거합니다.
- B도 매칭에서는 수관 polygon의 envelope 상자를 사용하므로 네 방법 모두 동일한 point-in-box 헝가리안 규칙입니다.
- `method_comparison.csv` 한 표에 네 방법×전체/GRSM/ORNL F1·각 단계 Bias·탄소 평균·비율·r·RMSE/MAE를 저장합니다.
- 23개 공통 조사구 탄소 RMSE 최소 방법을 이후 단계 기본값으로 연결합니다. 같은 표본에서 선택한 점수이며 새 지역 독립 검증은 아닙니다.
- 네 요청 조사구의 `*_methods_side_by_side.png`, `A_precision_recall.png`, fold 분할·선택 기록을 저장합니다.
- 실제 기본 방법은 Colab 실행 후 `selected_detection_method.json`과 report.md 맨 위에 기록됩니다.

ORNL 누락 원인: 공식 products API에는 RGB·CHM 모두 2021-06이 있지만 **PROVISIONAL**로 분류되어 있었습니다.
이제 해당 자료를 명시적으로 포함합니다. ORNL 다운로드 좌표는 UTM 16N, 작업 CRS는 17N을 사용합니다.

## 사용자가 할 일

1. Google Colab의 파일 → 노트북 업로드에서 ipynb 파일을 선택합니다.
2. 런타임 → 런타임 유형 변경 → T4 GPU를 선택합니다.
3. 첫 패키지 설치 셀을 실행합니다. 재시작 안내가 나오면 재시작합니다.
4. NEON 계정의 API Tokens에서 토큰을 발급해 Colab 보안 비밀 `NEON_TOKEN`에 등록하거나 숨김 입력창에 입력합니다.
5. 설정 셀을 실행하고 Drive 접근을 허용합니다. 기본은 GRSM·ORNL, 2021년 실측만 사용합니다.
6. 이번 비교는 `max_plots_per_site=None`을 유지하세요. 같은 23개 조사구를 모두 사용해야 합니다.
7. 위에서 아래로 실행합니다. 같은 파일·설정·패키지 버전의 완료 다운로드와 탐지는 재사용합니다.
8. `outputs/실행ID/audit.csv`와 조사구 overlay 그림부터 확인합니다. 다음으로 오차표·조사구 표·report.md를 봅니다.
9. 2021년 실측 조사구가 부족하면 `allow_adjacent_years=True`로 바꿉니다. 보완 연도는 결과에 따로 표시됩니다.
10. 영상 연도 대체는 별도 설정입니다. `imagery_allow_adjacent_years`는 ORNL만 True이며,
    2021 공통 자료가 없을 때에만 2022→2020을 검토합니다. 현재 제공 상태에서는 ORNL도 2021을 사용합니다.
11. 마지막 요약표를 확인한 뒤 ZIP 셀을 실행하면 `/content/neon_results.zip`을 다운로드합니다.
    Drive의 `outputs/neon_results.zip`에도 사본을 남깁니다.

모든 데이터와 결과의 Drive 루트는 `/content/drive/MyDrive/neon_grsm_ornl/`입니다.
예상 시간은 최초 1–수 시간이며 다운로드가 대부분을 차지할 수 있습니다.

## 셀 목록

| 절 | 내용 |
|---|---|
| 도입 | 실행 순서, 예상 시간, 분석 범위 |
| 01 | 필요한 패키지 설치 |
| 02 | Drive 마운트, 설정, 토큰 입력 |
| 03 | 공통 라이브러리, 다운로드/API 캐시, 버전 기록 |
| 04 | NEON 실측 다운로드와 연도 선택 |
| 05 | 기준점으로 줄기 좌표 계산, 실제 조사 소구획 경계 검산 |
| 06 | NEONForestAGB 참조값과 수종별 계수 로딩 |
| 07 | Jucker·Chojnacky 공식, 단위·방위각 검산 |
| 08 | API 제공 연월·release 조회, PROVISIONAL 포함, 영상 선택, 타일 다운로드·파일 수 로그 |
| 09 | 영상 색인, 조사구별 RGB/CHM 커버 여부 출력, 주변 잘라내기 |
| 10 | DeepForest 사전학습 모델 불러오기와 버전 호환 |
| 11 | 탐지, 수고·수관 지름·DBH·AGB·탄소 계산 |
| 12 | 보이는 나무 정의와 헝가리안 1:1 매칭 |
| 13 | 조사구 GroupKFold DBH 회귀, 수종 민감도 |
| 14 | 단계별 오차, 조사구 총량, 누락·미매칭 탐지 차이 분해 |
| 14A | A/B/C 함수 정의 (CHM 최대점·watershed·결합·중첩 GroupKFold) |
| 14B | 합성 CHM/정답 격리 점검, 15조합 탐색, 네 방법 평가, 나란히 그림, 기본값 선택 |
| 15 | 수고 ±10%, ±20%에 따른 조사구 탄소 민감도 |
| 16 | 조사구 RGB 그림과 1:1 산점도 |
| 17 | CSV, 진단, 한국어 report.md 저장 |
| 18 | 재사용 함수 neon_pipeline.py 내보내기 |
| 19 | 한국 자료에 적용할 때 바꿀 입력과 주의점 (설명만) |
| 20 | 사이트별 최종 요약표 (보고서 맨 위와 동일) |
| 21 | `/content/neon_results.zip` 저장·다운로드 |

## 중요한 해석

- 기본 모집단은 살아 있는 DBH≥10cm 나무입니다. 작은 나무의 별도 조사 면적을 섞지 않습니다.
- “전체 나무 기준” 탐지 지표도 이 모집단의 좌표가 확인된 나무에 한정됩니다.
- 조사 면적이 복원되지 않으면 해당 조사구를 제외하고 사유를 저장합니다. 빈 소구획을 추측하지 않습니다.
- 조사구 총량은 좌표 필터 전 `field_all`을 합칩니다. 좌표 누락은 총량을 NaN으로 만드는 조건이 아닙니다.
- `field_trees_all.csv`는 전체 참조 모집단, `field_trees.csv`는 매칭 가능한 좌표 있는 나무입니다.
- `AGB_ref_coordonly_t_ha`도 보존합니다. 실제 참조 AGB가 누락된 경우만 전체 총량을 NaN으로 둡니다.
- 분모는 `totalSampledAreaTrees`를 우선하고 `polygon_area_m2`·`area_difference_m2`로 면적 차이를 기록합니다.
- plot_AGB/plot_carbon 지표에 평균 추정·평균 참조·추정/참조 비율을 추가합니다.
- `plot_carbon_2021.png`는 x=전체 실측 기반 탄소, y=추정 탄소이며 사이트별 색을 사용합니다.
- `site_summary.csv`와 report.md 맨 위에 F1·각 단계 Bias·전체 실측 수·평균 탄소·비율을 표시합니다.
- Chojnacky에 “mixed hardwood”라는 단일 그룹은 없습니다. 일반 활엽수 대체 계수를 기본으로 표시하고 실측 수종별 식과 비교합니다.
- NEONForestAGB는 실측 DBH로 계산한 참조 추정값이며, 나무 무게를 직접 잰 값이 아닙니다.
- QUDO 등 충돌 수종 계수는 경고 후 제외하며 수종 기지 민감도에서 결측으로 표시합니다.
- gzip 등 Content-Encoding 응답의 크기는 전송 Content-Length와 비교하지 않습니다.
- HF 모델 캐시는 `/content/hf_cache`이며 Drive에 두지 않습니다.
- 재사용 모듈은 한국 자료가 있는 Colab에서 `configure(config, output_directory, projected_crs)` 후 사용하세요.

문법·구조 검증은 `VALIDATION.md`, 기존 GRSM 집계 검산은 `../outputs/revision_check/README.md`를 참고하세요.
기존 outputs/ 원본은 유지했습니다. 수정본의 ORNL 다운로드·GPU 탐지는 Colab에서 재실행해야 합니다.
나무 찾기 개선 추가 추론은 3 patch×23 조사구이며, threshold 5개는 동일 raw 예측을 필터링해 재사용합니다.
추가 처리 시간은 GPU·Drive 속도에 따라 수십 분에서 수 시간이 걸릴 수 있습니다.
