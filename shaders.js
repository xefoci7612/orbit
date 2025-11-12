'use strict';

// Inject shader code before compile and preprocessing, shader is still
// a list of #include directives.
//
// Original sources from:
// https://github.com/mrdoob/three.js/blob/dev/src/renderers/shaders/ShaderChunk/

export function addEarthShaders(shader, lightDirection, earthCloudMap, cloudsOffset) {

  const uniforms = {
    ulightDirection: { value: lightDirection },
    tClouds: { value: earthCloudMap },
    uv_xOffset: { value: cloudsOffset },
  };
  const atTop = `
    uniform vec3 ulightDirection;
    uniform sampler2D tClouds;
    uniform float uv_xOffset;
  `;
  const anchor = '#include <emissivemap_fragment>';
  const earthShader = `
    // --- Atmospheric fresnel effect ---

    vec3 viewDirection = normalize(vViewPosition);
    vec3 blueColor = vec3( 0.3, 0.6, 1.0 );

    // Fresnel effect color is due to Rayleigh Scattering and increases
    // with atmospheric path length
    float intensity = 1.1 - clamp(dot( vNormal, viewDirection ), 0.1, 1.0);
    vec3 cameraFresnel = blueColor * pow(intensity, 5.0);

    // Intensity fades out on the dark side
    float sunIntensity = smoothstep(0.0, 0.3, dot(vNormal, ulightDirection));

    vec3 fresnelEffect = cameraFresnel * sunIntensity;
    diffuseColor.rgb += fresnelEffect;

    // --- Cloud shadows on Earth ---

    // Our goal here is to use a “negative light map” approach to cast cloud shadows,
    // For any uv point X on earth map, we find the corresponding uv point on clouds
    // map that is directly above earth point X.
    //
    // Since the clouds spin eastward around earth (same increasing direction of uv.x),
    // we need the angular offset between cloud and earth, then, to get the correct cloud
    // texture point in this earth's fragment shader we minus Earth's uv.x by uv_xOffset.
    vec2 cloud2DPoint = vec2(vMapUv.x - uv_xOffset, vMapUv.y);
    float cloudIntensity = texture2D(tClouds, cloud2DPoint).r;

    // A more intense shadow is achieved by multiplying the surface color by
    // a smaller number. We also clamp the shadowValue to a minimum so it
    // doesn't get too dark
    diffuseColor.rgb *= max(1.0 - cloudIntensity, 0.1);
  `;
  Object.assign(shader.uniforms, uniforms);
  shader.fragmentShader = atTop + '\n' + shader.fragmentShader;
  shader.fragmentShader = shader.fragmentShader.replace(anchor, anchor + earthShader);
}

export function addCloudsShaders(shader, lightDirection) {

  const uniforms = {
    ulightDirection: { value: lightDirection },
  };
  const atTop = `
    uniform vec3 ulightDirection;
  `;
  const anchor1 = '#include <emissivemap_fragment>';
  const cloudShader = `

    #ifdef USE_EMISSIVEMAP

    	vec4 emissiveColor = texture2D( emissiveMap, vEmissiveMapUv );

    	#ifdef DECODE_VIDEO_TEXTURE_EMISSIVE

    		// use inline sRGB decode until browsers properly support SRGB8_ALPHA8 with video textures (#26516)

    		emissiveColor = sRGBTransferEOTF( emissiveColor );

    	#endif

      // --- Day/Night emissive logic ---
      //
      // NdotL is 1 on the sun side and 0 on the dark side, we normalize
      // vNormal due to interpolation from vertex to pixel shader
      float NdotL = dot(normalize(vNormal), ulightDirection);

      // smoothstep yields 1.0 on the day side and 0.0 on the night side
      float dayFactor = smoothstep(-0.2, 0.2, NdotL);
      emissiveColor.rgb *= dayFactor;

    	totalEmissiveRadiance *= emissiveColor.rgb;

    #endif
  `;
  Object.assign(shader.uniforms, uniforms);
  shader.fragmentShader = atTop + '\n' + shader.fragmentShader;
  shader.fragmentShader = shader.fragmentShader.replace(anchor1, cloudShader);

  const anchor2 = '#include <opaque_fragment>';
  const darkShader = `
    // Multiply outgoingLight (the combined ambient, diffuse, and specular light) by
    // our previously defined dayFactor to make the clouds black on the night side.
    outgoingLight.rgb *= dayFactor;
  `;
  shader.fragmentShader = shader.fragmentShader.replace(anchor2, darkShader + '\n' + anchor2);
}
