'use strict';

// Inject shader code before compile and preprocessing, shader is still
// a list of #include directives.
//
// Original sources from:
// https://github.com/mrdoob/three.js/blob/dev/src/renderers/shaders/ShaderChunk/

export function addFresnelShader(shader, lightDirection) {

  const uniforms = {
    ulightDirection: { value: lightDirection }
  };
  const atTop = `
    uniform vec3 ulightDirection;
  `;
  const anchor = '#include <emissivemap_fragment>';
  const fresnelShader = `
    // --- Atmospheric fresnel effect ---

    vec3 viewDirection = normalize(vViewPosition);
    vec3 blueColor = vec3( 0.3, 0.6, 1.0 );

    // Fresnel effect color is due to Rayleigh Scattering and increases
    // with atmospheric path length
    float intensity = 1.0 - clamp(dot( vNormal, viewDirection ), 0.1, 1.0);
    vec3 cameraFresnel = blueColor * pow(intensity, 5.0);

    // Intensity fades out on the dark side
    float sunIntensity = smoothstep(0.0, 0.3, dot(vNormal, ulightDirection));

    vec3 fresnelEffect = cameraFresnel * sunIntensity;
    diffuseColor.rgb += fresnelEffect;
  `;
  Object.assign(shader.uniforms, uniforms);
  shader.fragmentShader = atTop + '\n' + shader.fragmentShader;
  shader.fragmentShader = shader.fragmentShader.replace(anchor, anchor + fresnelShader);
}
