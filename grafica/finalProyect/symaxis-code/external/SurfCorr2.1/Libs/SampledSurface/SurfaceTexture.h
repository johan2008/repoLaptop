#ifndef __SURFACE_TEXTURE_H
#define __SURFACE_TEXTURE_H

#include "gaps.h"
#include "ParamParser.h"

class SurfaceTexture
{
	public:
		SurfaceTexture(const VKString & textureName, R2Image * image);
		virtual ~SurfaceTexture();
		virtual int InitializeTexture();
		virtual int LoadTexture();
		
		virtual void Draw(ParamParser & drawParams);
	protected:
		VKString m_TextureName;
		bool m_bInitializeTexture;
		R2Image * m_pImage;
		unsigned int m_TextureID;
		static unsigned int s_MaxTextureID;
};

#endif


