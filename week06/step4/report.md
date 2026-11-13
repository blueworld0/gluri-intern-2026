# Week06 Step 4 — CDSE 인증 및 다운로드 준비 보고서

점검일: 2026-11-13 (한국 시간)

준비 점검 당시 인증 도구와 저장 경로를 준비하고 공개 제품 메타데이터와 OIDC 설정만 조회했다. 아래 준비 점검 표는 당시 상태다. 이후 사용자가 사이트에서 직접 다운로드한 파일을 확인했으며, 후속 결과는 다음과 같다. Week03~05와 Step2~3 자료는 수정하지 않았다.

## 사용자 다운로드 후 확인 결과

2026-11-13 18:13 한국 시간에 ZIP을 읽기 전용으로 검사했다.

- 대상 제품명과 XML의 제품명이 일치한다.
- 파일 크기 974,110,105 bytes 및 MD5가 공식 카탈로그 기록과 일치한다.
- ZIP 전체 CRC 검사 통과, 필수 25개 파일 모두 포함, 파일 크기 일치.
- 원제품 XML에서 BOA_QUANTIFICATION_VALUE=10000, 밴드별 BOA_ADD_OFFSET=-1000을 확인했다.
- 이 원제품의 반사율은 NoData를 먼저 제외한 뒤 `(DN - 1000) / 10000`으로 계산한다. 기존 Analytical TIFF의 `/65535` 타당성이 입증됐다는 뜻은 아니다.
- ZIP 압축 해제와 영상 픽셀 처리는 아직 수행하지 않았다.

상세 기록: [다운로드 검증 결과](checks/20261113T091309Z_download_cddb9b9a/download_verification.json).

## 준비 상태

| 항목 | 확인 결과 |
| --- | --- |
| 인증 도구 | 준비 완료. 로컬 PowerShell의 숨김 입력으로 인증 확인 가능 |
| 실제 계정 인증 | 미실행 |
| 공개 제품 조회 | HTTP 200, UUID와 SAFE 제품명 일치 |
| 제품 상태 | Online=true |
| 공식 인증 설정 | 토큰 엔드포인트 일치, password grant 지원 확인 |
| 계정의 다운로드 권한 | 미확인. 인증된 다운로드 요청을 보내지 않음 |
| 실제 다운로드 | 시작하지 않음. 준비 도구에는 영상 다운로드 기능이 없음 |
| 기존 파일 충돌 | 대상 ZIP·SAFE 폴더·부분 다운로드 파일 모두 없음 |
| 저장 경로 | 쓰기 가능. downloads 폴더는 Git에서 제외 |

## 대상 제품과 용량

```text
S2C_MSIL2A_20250428T020711_N0511_R103_T52SDD_20250428T054014.SAFE
UUID: 6c4173d8-59e8-41f8-b339-c49620900f8f
```

| 구분 | 용량 |
| --- | --- |
| 전체 제품의 카탈로그 ContentLength | 974,110,105 bytes, 약 974.11MB |
| Step3에서 선별한 25개 파일 합계 | 485,041,677 bytes, 약 485.04MB |
| D: 드라이브 여유 공간 | 점검 당시 약 845.51GiB |
| 준비 단계 최소 여유 공간 기준 | 5GiB |

선별 파일 합계와 전체 SAFE 제품 용량은 다르다. 실제 ZIP 전송 크기와 압축 해제 용량은 다운로드 응답 및 압축 파일 목록을 확인한 뒤 확정한다. 5GiB는 ZIP·부분 파일·압축 해제 공간을 고려한 준비 기준이다.

저장 경로:

```text
D:\Coding\gluri-intern-2026\week06\step4\downloads\
```

향후 다운로드 시 같은 이름의 ZIP·SAFE 폴더·부분 파일이 있으면 덮어쓰지 않고 중단해야 한다. 이번에는 다운로드 파일을 생성하지 않았다.

## 로컬 환경과 공식 인증 방식

- Windows 11, PowerShell 5.1.26100.9818 / Desktop
- Python 3.12.2: `week03/.venv/Scripts/python.exe`
- OpenSSL 3.0.13, 기본 TLS 인증서 검증 사용
- Python 표준 라이브러리만 사용하므로 추가 설치 불필요
- 기존 가상환경 패키지는 변경하지 않음

[CDSE 공식 토큰 문서](https://documentation.dataspace.copernicus.eu/APIs/Token.html)에 따라 공식 HTTPS 엔드포인트에 `client_id=cdse-public`, `grant_type=password`, 사용자명과 비밀번호를 요청 본문으로 전송한다. 2FA 사용 계정은 현재 인증 코드인 `totp`도 보낸다. 토큰 유효기간은 응답의 `expires_in`으로 확인한다.

```text
https://identity.dataspace.copernicus.eu/auth/realms/CDSE/protocol/openid-connect/token
```

[공식 OData 다운로드 문서](https://documentation.dataspace.copernicus.eu/APIs/OData.html#product-download)에 따르면 실제 제품 다운로드에는 액세스 토큰을 담은 `Authorization: Bearer` 헤더가 필요하다. 대상은 `Products(UUID)/$value` 주소다. 이번에는 이 주소를 호출하지 않았다. 이 방식에는 별도의 S3 키나 Sentinel Hub client secret이 필요하지 않다.

## 인증정보 처리

준비한 도구는 사용자명·비밀번호·2FA 코드를 로컬 숨김 입력으로 받는다. 인증정보를 명령 인수, 환경변수, 설정 파일 또는 `.env`에 넣지 않는다. 비밀번호·토큰·인증 응답 본문을 파일이나 화면에 출력하지 않는다.

비대화형 실행과 숨김 입력 불가 상황에서는 중단한다. TLS 인증서 검증을 사용하고 인증 요청의 리다이렉트를 차단한다. 오류 발생 시 응답 본문이나 인증정보가 traceback으로 출력되지 않도록 처리했다.

인증 성공 여부와 만료시간만 표시한 뒤 종료한다. 토큰은 프로세스 메모리에서만 처리하므로 향후 다운로드 과정에서는 새로 인증해야 한다. Python 문자열 메모리의 완전한 덮어쓰기 삭제를 보장한다는 의미는 아니다.

## 사용자가 다음에 할 작업

본인 PC의 일반 PowerShell에서 다음 명령을 실행한다. 비밀번호나 토큰을 대화창에 입력하거나 명령 뒤에 붙이지 않는다.

```powershell
Set-Location 'D:\Coding\gluri-intern-2026'
& .\week03\.venv\Scripts\python.exe -B .\week06\step4\prepare_cdse.py --check-auth
```

1. 첫 번째 숨김 프롬프트에 CDSE 사용자명 또는 이메일을 입력한다.
2. 두 번째 숨김 프롬프트에 비밀번호를 입력한다.
3. 2FA를 사용하면 인증 앱의 현재 코드를 입력한다. 사용하지 않으면 Enter를 누른다.
4. `Authentication succeeded` 문구를 확인한다.

이 명령은 인증 확인만 수행한다. 다운로드를 시작하거나 토큰을 파일에 저장하지 않는다. 입력을 별도로 수집하거나 녹화하는 터미널 세션은 사용하지 않는다. 결과 공유가 필요하면 성공·실패 또는 HTTP 상태만 전달한다.

로컬 준비 상태만 점검하려면 다음 명령을 사용한다.

```powershell
& .\week03\.venv\Scripts\python.exe -B .\week06\step4\prepare_cdse.py
```

공개 메타데이터까지 다시 확인하려면 다음 명령을 사용한다.

```powershell
& .\week03\.venv\Scripts\python.exe -B .\week06\step4\prepare_cdse.py --public-preflight
```

준비 점검 결과는 `checks/` 아래의 고유 디렉터리에 저장되어 이전 보고서를 덮어쓰지 않는다. 인증 확인 모드는 인증 응답을 저장하지 않는다.

## 검증 및 산출물

합성 입력만 사용하는 오프라인 보안 검사 5개를 통과했다. 실제 계정 인증 결과는 아니다.

1. 성공 출력에 사용자명·비밀번호·토큰이 포함되지 않음
2. HTTP 오류 응답의 본문과 헤더가 노출되지 않음
3. 비대화형 실행에서 입력 및 네트워크 요청 전에 인증을 차단함
4. 인증 요청의 리다이렉트를 따르지 않음
5. 공개 점검 함수에서 토큰 및 다운로드 주소 호출을 차단함

| 파일 | 내용 |
| --- | --- |
| [prepare_cdse.py](prepare_cdse.py) | 준비 점검 및 로컬 숨김 인증 확인 도구 |
| [product.json](product.json) | 비밀정보 없는 대상 제품과 공식 URL 설정 |
| [test_auth_safety.py](test_auth_safety.py) | 합성 입력을 사용하는 오프라인 보안 검사 |
| [validation.json](validation.json) | 검사 결과와 실제 인증·다운로드 미실행 상태 |
| [preflight.json](checks/20261113T045634Z_4eb39837/preflight.json) | 실제 로컬 및 공개 메타데이터 점검 결과 |
| [product_metadata.json](checks/20261113T045634Z_4eb39837/product_metadata.json) | 공개 제품 메타데이터 원응답 |
| [oidc_public_metadata.json](checks/20261113T045634Z_4eb39837/oidc_public_metadata.json) | 공개 OIDC 설정 원응답 |

준비 점검에서는 실제 계정 인증을 수행하지 않았다. 이후 **사용자가 직접 받은 ZIP의 제품 식별·용량·MD5·CRC·필수 파일 검증까지 완료**했다. API 인증 확인 없이도 후속 로컬 분석을 진행할 수 있다.
