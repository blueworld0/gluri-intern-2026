// (선택) Google Earth Engine 코드 에디터용: ESA CCI Biomass v6.0 (2007, 2010, 2015–2022)
// 파이썬 파이프라인은 v4(2017–2020)를 쓰는데, v6 가 더 최신·연도도 많습니다.
// 사용법: code.earthengine.google.com → 새 스크립트에 붙여넣기 →
//        왼쪽 Assets 에 경계(geojson/shp)를 올리고 아래 BOUNDARY 경로를 바꾼 뒤 Run.
var BOUNDARY = ee.FeatureCollection('projects/YOUR_PROJECT/assets/geumjeongsan_proxy');
var R = 0.43, CF = 0.5;
var tree = ee.ImageCollection('ESA/WorldCover/v100').first().eq(10);   // 2020 나무
var cci = ee.ImageCollection('ESA/CCI/Above_Ground_Biomass/V6_0');
var px_ha = ee.Image.pixelArea().divide(1e4);

var rows = cci.map(function (img) {
  var agb_t = img.select('agb').multiply(px_ha).updateMask(tree);       // 픽셀별 AGB 톤
  var s = agb_t.reduceRegion({reducer: ee.Reducer.sum(), geometry: BOUNDARY.geometry(),
                              scale: 10, maxPixels: 1e10});
  var agb = ee.Number(s.get('agb'));
  return ee.Feature(null, {year: img.date().get('year'), agb_t: agb,
                           carbon_tC: agb.multiply(1 + R).multiply(CF)});
});
print(rows);
Export.table.toDrive({collection: rows, description: 'geumjeongsan_cci_v6', fileFormat: 'CSV'});
Map.centerObject(BOUNDARY, 12);
Map.addLayer(cci.filter(ee.Filter.calendarRange(2022, 2022, 'year')).first().select('agb'),
             {min: 0, max: 200, palette: ['white', 'green']}, 'AGB 2022');
Map.addLayer(BOUNDARY, {color: 'black'}, '경계');
