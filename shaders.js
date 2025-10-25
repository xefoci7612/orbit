'use strict';

import * as THREE from 'three';

export class Shaders {

  // Blocks emissive map during day time and add fresnel effect
  //
  // https://github.com/mrdoob/three.js/blob/dev/src/renderers/shaders/ShaderChunk/emissivemap_fragment.glsl.js
  emissivemap(normalizedSunRay) {
    return {
      uniforms: {
        uSunRay: { value: normalizedSunRay },
      },
      fragmentShader: {
        atTop: `
          uniform vec3 uSunRay;
        `,
        // Anchors work before code is preprocessed, so it's a long list
        // of #include, we have to live with that
        anchor: '#include <emissivemap_fragment>',
        injection: `
        #ifdef USE_EMISSIVEMAP

        	vec4 emissiveColor = texture2D( emissiveMap, vEmissiveMapUv );

        	#ifdef DECODE_VIDEO_TEXTURE_EMISSIVE

        		// use inline sRGB decode until browsers properly support SRGB8_ALPHA8 with video textures (#26516)

        		emissiveColor = sRGBTransferEOTF( emissiveColor );

        	#endif

          // --- Day/Night emissive logic ---
          //
          // uSunRay is the normalized vector from the Sun to the Earth.
          // For standard lighting, we need the direction from the surface to the light.
          // This is the inverse of the sun ray's direction.
          float NdotL = dot(vNormal, -uSunRay);

          // smoothstep yields 1.0 on the day side and 0.0 on the night side
          float nightFactor = 1.0 - smoothstep(-0.5, 0.0, NdotL);
          emissiveColor.rgb *= nightFactor;

        	totalEmissiveRadiance *= emissiveColor.rgb;

          // --- Atmospheric fresnel effect ---
          //
          // adding a small amount of atmospheric fresnel effect to make it more realistic
          // fine tune the first constant below for stronger or weaker effect
          vec3 blueColor = vec3( 0.3, 0.6, 1.0 );
          vec3 uToCamera = vec3( 0.0, 0.0, 1.0 );
          float intensity = 1.0 - clamp(dot( vNormal, uToCamera ), 0.2, 1.0);
          vec3 fresnel = blueColor * pow(intensity, 5.0);

          diffuseColor.rgb += fresnel;

        #endif
        `
      }
    };
  };

  injectGLSL(shader, chunks) {

    let newSource = shader.fragmentShader;

    chunks.forEach((chunk, index) => {
      if (chunk.uniforms) {
        Object.assign(shader.uniforms, chunk.uniforms);
      }

      const sh = chunk.fragmentShader;

      if (sh.atTop)
        newSource = sh.atTop + '\n' + newSource;

      if (!newSource.includes(sh.anchor))
        console.warn(`[Chunk ${index}]: ANCHOR NOT FOUND. Injection failed for:\n"${sh.anchor}"`);

      newSource = newSource.replace(sh.anchor, sh.injection);
    });

    shader.fragmentShader = newSource;
  }
};


// Creates a complete ShaderMaterial configuration for a physically-based
// atmospheric scattering shell.
export function atmosphereShaderConfig(normalizedSunRay) {

  const vertexShader = `
    // We receive the vertex normal from the geometry
    varying vec3 vNormal;

    void main() {
      // Pass the normal to the fragment shader, transformed into View Space
      vNormal = normalize( normalMatrix * normal );

      // Calculate the final position of the vertex
      gl_Position = projectionMatrix * modelViewMatrix * vec4( position, 1.0 );
    }
  `;

  const fragmentShader = `
    // We receive the View Space normal from the vertex shader
    varying vec3 vNormal;

    // We will pass the sun's direction in View Space as a uniform
    uniform vec3 uSunDirectionView;

    // Blue gradient for atmosphere height dependent color
    uniform sampler2D uGradientMap;

    void main() {
      // 1. FRESNEL (Limb) EFFECT
      // The dot product gives us a value from -1 (center) to 0 (edge).
      // We use smoothstep to create a nice falloff from 0.0 to 1.0 at the limb.
      float fresnel = smoothstep( 0.7, -0.3, vNormal.z );

      // 2. SUN SCATTERING EFFECT
      // We get the dot product of the normal with the sun's direction.
      float sunDot = dot( vNormal, -uSunDirectionView );

      // We only want the glow on the sunlit side, so we clamp negative values.
      float sunlitFactor = max( sunDot, 0.0 );

      // pow() gives us a bright "hotspot" where the sun is hitting directly,
      // and it fades out nicely.
      sunlitFactor = pow( sunlitFactor, 3.0 );

      // The final brightness is the sum of the general limb glow and the sun hotspot.
      float intensity = fresnel * 0.4 + sunlitFactor * 0.8;

      // The atmosphere gets whiter/brighter where the sun hits it most directly.
      vec3 baseColor = vec3(0.3, 0.6, 1.0); // Base sky blue
      vec3 hazeColor = vec3(0.7, 0.8, 1.2); // Still bright, but distinctly blue
      vec3 atmosphereColor = mix(baseColor, hazeColor, pow(sunlitFactor, 2.0));

      // The final color is the atmosphere color multiplied by our calculated intensity.
      // We set the alpha to the intensity to make it nicely transparent.
      gl_FragColor = vec4( atmosphereColor * intensity, intensity );
    }
  `;

  return {
    vertexShader: vertexShader,
    fragmentShader: fragmentShader,
    blending: THREE.AdditiveBlending,
    transparent: true,
    depthWrite: false,
    side: THREE.BackSide,
    uniforms: {
      uSunDirectionView: { value: normalizedSunRay }
  }
  };
};
