//VERSION=3
// DRAFT ONLY: not executed. Requires explicit authenticated Process API use.
function setup() {
  return {
    input: [
      {
        bands: ["B02", "B03", "B04", "B08", "B11", "B12", "SCL", "dataMask"],
        units: ["REFLECTANCE", "REFLECTANCE", "REFLECTANCE", "REFLECTANCE", "REFLECTANCE", "REFLECTANCE", "DN", "DN"]
      }
    ],
    output: [
      { id: "reflectance", bands: 6, sampleType: "FLOAT32" },
      { id: "quality", bands: 2, sampleType: "UINT8" }
    ],
    mosaicking: "TILE"
  };
}
function evaluatePixel(samples) {
  // Example only: single-scene validation requests must first verify the scenes.
  let s = samples[0];
  if (!s) return {reflectance: [0,0,0,0,0,0], quality: [0,0]};
  return {reflectance: [s.B02,s.B03,s.B04,s.B08,s.B11,s.B12], quality: [s.SCL,s.dataMask]};
}
function updateOutputMetadata(scenes, inputMetadata, outputMetadata) {
  outputMetadata.userData = { scenes: scenes, inputMetadata: inputMetadata };
}
