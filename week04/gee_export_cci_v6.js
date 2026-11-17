// week04 · ESA CCI Biomass v6.0 금정산 주변 잘라 받기 (Google Earth Engine 코드 에디터용)
// 사용법: code.earthengine.google.com → New Script → 전체 붙여넣기 → Run
//        → 오른쪽 [Tasks] 탭에서 연도별 RUN 클릭 → 구글 드라이브 'geumjeongsan_cci_v6' 폴더에 저장됨
// 결과: 연도마다 GeoTIFF 1개 (밴드 1 = agb, 밴드 2 = agb_sd, 단위 t/ha, 100 m 원래 격자 그대로)

var BBOX = ee.Geometry.Rectangle([128.98, 35.12, 129.14, 35.36]);   // week03 config.py 와 같은 범위
var YEARS = [2015, 2016, 2017, 2018, 2019, 2020, 2021, 2022];

var cci = ee.ImageCollection('ESA/CCI/Above_Ground_Biomass/V6_0');
print('v6 이미지 목록(연도 확인용)', cci.aggregate_array('system:index'));
print('밴드 이름', cci.first().bandNames());

var proj = cci.first().projection().getInfo();                      // 원래 격자를 그대로 유지

YEARS.forEach(function (y) {
  var img = cci.filterDate(y + '-01-01', (y + 1) + '-01-01').first()
               .select(['agb', 'agb_sd']).toFloat();
  Export.image.toDrive({
    image: img,
    description: 'cci_v6_' + y,
    folder: 'geumjeongsan_cci_v6',
    fileNamePrefix: 'cci_v6_' + y,
    region: BBOX,
    crs: proj.crs,
    crsTransform: proj.transform,
    maxPixels: 1e9
  });
});

// 눈으로 확인
Map.centerObject(BBOX, 12);
Map.addLayer(cci.filterDate('2022-01-01', '2023-01-01').first().select('agb'),
             {min: 0, max: 250, palette: ['white', 'yellow', 'green', 'darkgreen']}, 'AGB 2022');
Map.addLayer(ee.FeatureCollection([ee.Feature(BBOX)]).style({color: 'red', fillColor: '00000000'}), {}, '범위');
