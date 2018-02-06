#include "SurfaceTexture.h"
#include "PipelineGeneral.h"

unsigned int SurfaceTexture::s_MaxTextureID = 0;

SurfaceTexture::SurfaceTexture(const VKString & texturename, R2Image * image)
{
	m_TextureName = texturename;
	assert(image!=NULL);
	m_pImage = image;
	m_bInitializeTexture = true;
}

SurfaceTexture::~SurfaceTexture()
{
}

int SurfaceTexture::InitializeTexture()
{
	m_TextureID = s_MaxTextureID++;
	glBindTexture (GL_TEXTURE_2D, m_TextureID);
    glPixelStorei (GL_UNPACK_ALIGNMENT, 0);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexEnvf (GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
	if (m_pImage->NComponents()==1)
		glTexImage2D (GL_TEXTURE_2D, 0, 1, m_pImage->Width(), m_pImage->Height(), 
					  0, GL_LUMINANCE, GL_UNSIGNED_BYTE, m_pImage->Pixels());
	else if (m_pImage->NComponents()==3)
		glTexImage2D (GL_TEXTURE_2D, 0, GL_RGB, m_pImage->Width(), m_pImage->Height(), 
					  0, GL_RGB, GL_UNSIGNED_BYTE, m_pImage->Pixels());
	else
		assert(false);
	                    
	return m_TextureID;
}

int SurfaceTexture::LoadTexture()
{
//	std::cout<<"Loading texture"<<std::endl;
	if (m_bInitializeTexture)
	{
		InitializeTexture();
		m_bInitializeTexture = false;
	}
	glBindTexture( GL_TEXTURE_2D, m_TextureID );
	return m_TextureID;
}

void SurfaceTexture::Draw(ParamParser & drawParams)
{
	bool valid;
	glEnable(GL_TEXTURE_2D);
	glDisable(GL_LIGHTING);
	glColor3d(1., 1., 1.);
	
	LoadTexture();
	glBegin( GL_QUADS );
	glTexCoord2d(0.0,0.0); glVertex3d(0.0,0.0, -1.);
	glTexCoord2d(1.0,0.0); glVertex3d(1.0,0.0, -1.);
	glTexCoord2d(1.0,1.0); glVertex3d(1.0,1.0, -1.);
	glTexCoord2d(0.0,1.0); glVertex3d(0.0,1.0, -1.);
	glEnd();
	glDisable(GL_TEXTURE_2D);
	
	static GLfloat material[4];	
	material[3] = 1.;
	
	// Draw correspondences
	bool interactCorrs = drawParams.GetStrValue("RendererDefault", "InteractionMode", valid)=="ManualCorrespondences";
	
	if (interactCorrs)
	{
		SurfaceSampleSet * sampleSet = Surface2DPlane::m_pInfinitePlanePseudosurface->GetSampleSet(m_TextureName+AnalysisWindow::s_InteractiveSampleSetName);

		if (sampleSet!=NULL)
		{
			glEnable(GL_LIGHTING);
			double radius = .02;
			for (int i=0; i<sampleSet->NumSamples(); i++)
			{
				AnalysisWindow::s_pMainWindow->MapIntToColor(i, material);
				glMaterialfv(GL_FRONT, GL_AMBIENT_AND_DIFFUSE, material); 
				R3Point p = sampleSet->GetSample(i).GetPosition();
				p.Z(-.9);
				R3Sphere(p, radius).Draw();
			}				
		}
	}
}
