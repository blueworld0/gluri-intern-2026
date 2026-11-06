# Week05 Step3 — Sentinel-2 TIFF 원본 출력 설정·반사율 교차 검증

## 결론

**2025 TIFF의 `/65535` 변환은 C: 근거 부족으로 판단 불가**다. 저장소에는 2025 TIFF 생성 요청·evalscript·독립 FLOAT32 밴드 반사율이 없고, TIFF 자체에도 반사율 계수나 제품 ID가 없다. 현재 계산이 재현된다는 사실은 물리적 정확성을 입증하지 않는다. 실제 동시점·동일 픽셀 교차 검증은 수행할 수 없어 수치 비교 CSV·그래프를 만들지 않았다.

## 확인된 사실과 미확인 설정

- **확인:** 2025 TIFF 파일명은 Sentinel-2 L2A, 2025-04-28 요청 구간, B02/B03/B04/B08/B11/B12를 나타낸다. 첫 밴드는 `uint16`, 둘째는 0/1 마스크다. Week03 코드가 `1/65535`, offset 0을 **잠정 적용**했다(`week03/analyze_sentinel_2025.py:3-4,36,276`; `week05/step2/audit_metadata.json`).
- **출처 기록:** 2026 TIFF에는 Copernicus Browser Analytical 16-bit Raw Bands 사용을 확인한 설정 파일이 있다. 2025 TIFF는 유사 형식으로 취급했으나 2025 자체의 다운로드 요청 기록은 현재 저장소에서 미확인이다(`week03/reflectance_config.json:2`, `week03/README_SENTINEL_2025.md:13-16`).
- **미확인:** 2025 TIFF의 실제 생성일, 제품 ID, 정확한 촬영 시각, API 요청 본문, evalscript, `units`, `sampleType`, `harmonizeValues`, 반환식·출력 격자 설정, 추가 배율·offset·클리핑 여부. SCL JP2 이름의 촬영 시각과 TIFF가 동일 장면인지도 독립 검증 불가다(`week03/outputs/sentinel_2025/scl_validation/scl_quality_report.json`의 `tiff_association`).
- **독립 기준 부재:** 저장소의 FLOAT32 파일은 이미 잠정 변환으로 만든 지수 결과물이다. 2026 영상은 다른 시점이다. 어느 쪽도 2025 밴드 반사율의 직접 검증 기준이 아니다.

## 발견된 문제와 필요한 자료

가장 큰 문제는 **출력 계보의 단절**이다. TIFF의 `scale=1`, `offset=0` 및 빈 단위 태그로는 `REFLECTANCE ×65535` 여부를 알 수 없다. 원제품 L2A DN의 공식 offset·quantification 규칙과 Browser TIFF의 UINT16 출력 규칙도 동일시할 수 없다. 동일 장면에서 실제 설정을 확인하려면 (1) 2025 Browser 다운로드·API 요청과 evalscript, (2) 동일 제품의 SAFE/제품 XML 및 TIFF–SCL 연결 정보, (3) 같은 장면·격자의 검증된 FLOAT32 밴드 반사율이 필요하다. 자세한 근거와 비교 절차는 [출력 설정 추적](source_trace.md), [교차 검증 조사](reflectance_crosscheck.md), [변환 판정](conversion_assessment.md)에 기록했다.

## Step4 모델 개선 실험 진행 여부

**현재 수치로 개선 효과를 주장하는 실험은 보류하는 것이 타당하다.** 입력의 방사 값이 미검증이며, 특히 EVI가 공통 배율에도 민감하다. Step4에는 기존 점수·분할의 재현 확인이나 공간 검증 설계 등 **데이터를 바꾸지 않는 준비 작업**은 가능하다. 새 모델의 성능을 비교·해석하려면 먼저 출력 변환을 검증하고, 필요 시 별도 버전의 지수를 만들어 동일 평가 설계로 비교해야 한다. 이번 Step3에는 모델 재학습과 기존 데이터 수정을 수행하지 않았다.

## 질문에 대한 답

1. **기존 TIFF 생성 설정:** 2025 파일명에서 L2A·밴드·요청 날짜, 파일에서 `uint16`·마스크·격자만 확인된다. 2025 출력 요청의 세부 설정은 **현재 저장소에서 미확인**이다. 2026 파일의 Browser 사용 기록을 2025의 확정 설정으로 간주할 수 없다.
2. **`/65535` 타당성:** **판단 불가(C)**. 타당함(A)이나 오류(B)를 입증하는 실제 생성식·독립 반사율이 없다.
3. **추가 데이터:** **필요하다.** 2025 요청·evalscript, 제품 메타데이터, 동일 장면 FLOAT32 반사율을 우선 확보해야 한다.
4. **현재 데이터로 모델 개선 실험:** 탐색적 준비는 가능하지만, **신뢰할 수 있는 성능 개선 결론을 내리는 실험은 보류**한다. 방사 변환 검증 후 재평가한다.

## API 사용·비용

외부 데이터 API는 호출하지 않았다. [CDSE의 처리 단위 문서](https://documentation.dataspace.copernicus.eu/APIs/SentinelHub/Overview/ProcessingUnit.html)에 따르면 성공한 Sentinel Hub Process API 요청은 PU를 차감하고 FLOAT32 출력은 추가 계수가 있다. 향후 검증 요청은 계정 설정과 예상 사용량 확인이 선행돼야 한다.
