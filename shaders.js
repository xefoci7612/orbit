'use strict';

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
          float nightFactor = 1.0 - smoothstep(-0.3, 0.0, NdotL);
          emissiveColor.rgb *= nightFactor;

        	totalEmissiveRadiance *= emissiveColor.rgb;

          // --- Atmospheric fresnel effect ---
          //
          // adding a small amount of atmospheric fresnel effect to make it more realistic
          // fine tune the first constant below for stronger or weaker effect
          vec3 blueColor = vec3( 0.3, 0.6, 1.0 );
          vec3 uToCamera = vec3( 0.0, 0.0, 1.0 );
          float intensity = 1.0 - clamp(dot( vNormal, uToCamera ), 0.2, 1.0);
          vec3 atmosphere = blueColor * pow(intensity, 5.0);

          diffuseColor.rgb += atmosphere;

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
