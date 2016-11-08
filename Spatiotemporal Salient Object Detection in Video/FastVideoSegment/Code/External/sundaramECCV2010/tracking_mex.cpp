#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include "CTensor.h"
#include "CTensorMatlab.h"
#include "CFilter.h"

#include <matrix.h>
#include <mex.h>

class CTrack {
public:
	CTrack() {mStopped = false; mLabel = -1;}
	std::vector<float> mx,my;						 // current position of the track
	int mox,moy;													// original starting point of the track
	int mLabel;													 // assignment to a region (ignored for tracking but makes the tracking files compatible to other tools I have)
	bool mStopped;												// tracking stopped due to occlusion, etc.
	int mSetupTime;											 // Time when track was created
};

CVector<std::string> mInput;
std::string mFilename;
std::string mInputDir;
std::string mResultDir;
std::vector<CTrack> mTracks;
int mStartFrame;
int mStep;
int mSequenceLength;
int mXSize,mYSize;
CMatrix<float> mColorCode;
std::ofstream mStatus;

// buildColorCode
void buildColorCode();
void computeCorners( CTensor<float> &aImage, CMatrix<float> &aCorners, float aRho );
void dt( CVector<float> &f, CVector<float> &d, int n );
void euclideanDistanceTransform( CMatrix<float> &aMatrix );
void writeTracks();

class CPoint
{
	public:
		CPoint() {}
		float x,y,frame;
};

class CSimpleTrack
{
	public:
		CSimpleTrack() {}
		int mLabel;
		CVector<CPoint> mPoints;
};

void mexFunction( int nlhs, mxArray *plhs[], int nrhs,
	const mxArray *prhs[] )
{
	const mxArray *imgs;
	const mxArray *fflow;
	const mxArray *bflow;
	
	unsigned int frames;
	
	imgs = prhs[ 0 ];
	fflow = prhs[ 1 ];
	bflow = prhs[ 2 ];

	frames = mxGetNumberOfElements( imgs );

	CVector<float>color = CVector<float>(3);
	
	// Set beginning and end of tracking
	mStartFrame = 0;
	mSequenceLength = frames;
	mStep = 2;

	clock_t start = clock();
	// Load first image
	CTensorMatlab<float>* aImage1 = new CTensorMatlab<float>;
	CTensorMatlab<float>* aImage2;

	aImage1->copyFromMxArrayImage( mxGetCell( imgs, 0 ) );

	mXSize = aImage1->xSize();
	mYSize = aImage1->ySize();

	mTracks.clear();
	CMatrix<float> aCorners;
	CMatrix<float> aCovered(mXSize,mYSize);
	int aSize = mXSize*mYSize;
	// Smooth first image (can be removed from optical flow computation then)
	NFilter::recursiveSmoothX( *aImage1, 0.8f );
	NFilter::recursiveSmoothY( *aImage1, 0.8f );
	// Tracking
	for (int t = 0; t < frames - 1; t++)
	{
		// Load next image
		aImage2 = new CTensorMatlab<float>;
		aImage1->copyFromMxArrayImage( mxGetCell( imgs, t + 1 ) );

		NFilter::recursiveSmoothX( *aImage2, 0.8f );
		NFilter::recursiveSmoothY( *aImage2, 0.8f );
		// Mark areas sufficiently covered by tracks
		aCovered = 1e20;
		if (t > 0)
		{
			for (unsigned int i = 0; i < mTracks.size(); i++)
			{
				if (mTracks[i].mStopped == false)
					aCovered( ( int )mTracks[ i ].mx.back(), ( int )mTracks[ i ].my.back() ) = 0.0f;
			}
			euclideanDistanceTransform(aCovered);
		}
		
		// Set up new tracking points in uncovered areas
		computeCorners( *aImage1, aCorners, 3.0f );
		float aCornerAvg = aCorners.avg();
		for (int ay = 4; ay < mYSize-4; ay+=mStep)
		{
			for (int ax = 4; ax < mXSize-4; ax+=mStep)
			{
				if (aCovered(ax,ay) < mStep*mStep)
					continue;
				float distToImageBnd = exp( -0.1 * NMath::min( NMath::min( NMath::min( ax, ay ), mXSize - ax ), mYSize - ay ) );
				if( aCorners( ax, ay ) < 1.0 * ( aCornerAvg * ( 0.1 + distToImageBnd ) ) )
					continue;
				if( aCorners( ax, ay ) < 1.0 * ( 1.0f + distToImageBnd ) )
					continue;
					
				mTracks.push_back( CTrack() );
				CTrack& newTrack = mTracks.back();
				newTrack.mox = ax;
				newTrack.moy = ay;
				newTrack.mLabel = -1;
				newTrack.mSetupTime = t;
			}
		}
		
		// Compute bidirectional LDOF or read from file when available
		CTensorMatlab<float> aForward, aBackward;

		aForward.copyFromMxArrayOpticalFlow( mxGetCell( fflow, t ) );
		aBackward.copyFromMxArrayOpticalFlow( mxGetCell( bflow, t ) );
		
		// Check consistency of forward flow via backward flow
		CMatrix<float> aUnreliable( mXSize, mYSize, 0 );
		CTensor<float> dx( mXSize, mYSize, 2 );
		CTensor<float> dy( mXSize, mYSize, 2 );
		CDerivative<float> aDev( 3 );
		NFilter::filter( aForward, dx, aDev, 1, 1 );
		NFilter::filter( aForward, dy, 1, aDev, 1 );
		CMatrix<float> aMotionEdge( mXSize, mYSize, 0 );
		
		for( int i = 0; i < aSize; i++)
		{
			aMotionEdge.data()[ i ] += dx.data()[ i ] * dx.data()[ i ];
			aMotionEdge.data()[ i ] += dx.data()[ aSize + i ] * dx.data()[ aSize + i ];
			aMotionEdge.data()[ i ] += dy.data()[ i ] * dy.data()[ i ];
			aMotionEdge.data()[ i ] += dy.data()[ aSize + i ] * dy.data()[ aSize + i ];
		}
		
		int x1, y1, x2, y2;
		float bx, by, alphaX, alphaY, a, b, u, v, cx, cy, u2, v2;
		for( int ay = 0; ay < aForward.ySize(); ay++ )
		{
			for( int ax = 0; ax < aForward.xSize(); ax++ )
			{
				bx = ax+aForward(ax,ay,0);
				by = ay+aForward(ax,ay,1);
				x1 = floor(bx);
				y1 = floor(by);
				x2 = x1+1;
				y2 = y1+1;
				if( x1 < 0 || x2 >= mXSize || y1 < 0 || y2 >= mYSize )
				{
					aUnreliable(ax,ay) = 1.0f;
					continue;
				}
				
				alphaX = bx-x1;
				alphaY = by-y1;
				a = (1.0-alphaX)*aBackward(x1,y1,0)+alphaX*aBackward(x2,y1,0);
				b = (1.0-alphaX)*aBackward(x1,y2,0)+alphaX*aBackward(x2,y2,0);
				u = (1.0-alphaY)*a+alphaY*b;
				a = (1.0-alphaX)*aBackward(x1,y1,1)+alphaX*aBackward(x2,y1,1);
				b = (1.0-alphaX)*aBackward(x1,y2,1)+alphaX*aBackward(x2,y2,1);
				v = (1.0-alphaY)*a+alphaY*b;
				cx = bx+u;
				cy = by+v;
				u2 = aForward( ax, ay, 0 );
				v2 = aForward( ax, ay, 1 );
				
				if( ( ( cx - ax ) * ( cx - ax ) + ( cy - ay ) * ( cy - ay ) ) >= 0.01 * ( u2 * u2 + v2 * v2 + u * u + v * v ) + 0.5f )
				{
					aUnreliable( ax, ay ) = 1.0f;
					continue;
				}
				if( aMotionEdge( ax, ay ) > 0.01 * ( u2 * u2 + v2 * v2 ) + 0.002f )
				{
					aUnreliable( ax, ay ) = 1.0f;
					continue;
				}
			}
		}
		
		CTensor<float> aShow( *aImage2 );
		for( unsigned int i = 0; i < mTracks.size(); i++ )
		{
			if( mTracks[i].mStopped )
				continue;
			
			float ax, ay, oldvar;
			if( mTracks[i].mSetupTime == t )
			{
				ax = mTracks[ i ].mox;
				ay = mTracks[ i ].moy;
			}
			else
			{
				ax = mTracks[ i ].mx.back();
				ay = mTracks[ i ].my.back();
			}
			
			int iax = lroundf(ax);
			int iay = lroundf(ay);
			if( aUnreliable( iax, iay ) > 0 )
				mTracks[i].mStopped = true;
			else
			{
				float bx = ax + aForward( iax, iay, 0 );
				float by = ay + aForward( iax, iay, 1 );
				int ibx = lroundf( bx );
				int iby = lroundf( by );
				if( ibx < 0 || iby < 0 || ibx >= mXSize || iby >= mYSize )
					mTracks[ i ].mStopped = true;
				else
				{
					mTracks[ i ].mx.push_back( bx );
					mTracks[ i ].my.push_back( by );
					mTracks[ i ].mLabel = 0;
					int t = mTracks[ i ].mx.size();
					int f = 1300 / frames;
					f=10;
					color.data()[0] = mColorCode( 0, f * t );
					color.data()[1] = mColorCode( 1, f * t );
					color.data()[2] = mColorCode( 2, f * t );

					// draw colored square
					{
						int x1 = ibx-1 < 0 ? 0 : ibx-1;
						int y1 = iby-1 < 0 ? 0 : iby-1;
						int x2 = ibx+1 >= mXSize ? mXSize-1 : ibx+1;
						int y2 = iby+1 >= mYSize ? mYSize-1 : iby+1;
						aShow.fillRect(color, x1,y1,x2,y2);
					}
					
				}
			}
		}

		// Prepare for next image
		delete aImage1;
		aImage1 = aImage2;
	}
	delete aImage1; 

	return;
}

// buildColorCode
void buildColorCode() {
	mColorCode.setSize(3,15000);
	for (int i = 0; i < 256; i++) {
		mColorCode(0,i) = 0;
		mColorCode(1,i) = i;
		mColorCode(2,i) = 255;
	}
	for (int i = 0; i < 256; i++) {
		mColorCode(0,i+256) = 0;
		mColorCode(1,i+256) = 255;
		mColorCode(2,i+256) = 255-i;
	}
	for (int i = 0; i < 256; i++) {
		mColorCode(0,i+512) = i;
		mColorCode(1,i+512) = 255;
		mColorCode(2,i+512) = 0;
	}
	for (int i = 0; i < 256; i++) {
		mColorCode(0,i+768) = 255;
		mColorCode(1,i+768) = 255-i;
		mColorCode(2,i+768) = 0;
	}
	for (int i = 0; i < 256; i++) {
		mColorCode(0,i+1024) = 255;
		mColorCode(1,i+1024) = 0;
		mColorCode(2,i+1024) = i;
	}
	for (int i = 1280; i < mColorCode.ySize(); i++) {
		mColorCode(0,i) = 255;
		mColorCode(1,i) = 0;
		mColorCode(2,i) = 255;
	}
}

// computeCorners --------------------------------------------------------------
void computeCorners(CTensor<float>& aImage, CMatrix<float>& aCorners, float aRho) {
	aCorners.setSize(aImage.xSize(),aImage.ySize());
	int aXSize = aImage.xSize();
	int aYSize = aImage.ySize();
	int aSize = aXSize*aYSize;
	// Compute gradient
	CTensor<float> dx(aXSize,aYSize,aImage.zSize());
	CTensor<float> dy(aXSize,aYSize,aImage.zSize());
	CDerivative<float> aDerivative(3);
	NFilter::filter(aImage,dx,aDerivative,1,1);
	NFilter::filter(aImage,dy,1,aDerivative,1);
	// Compute second moment matrix
	CMatrix<float> dxx(aXSize,aYSize,0);
	CMatrix<float> dyy(aXSize,aYSize,0);
	CMatrix<float> dxy(aXSize,aYSize,0);
	int i2 = 0;
	for (int k = 0; k < aImage.zSize(); k++)
		for (int i = 0; i < aSize; i++,i2++) {
			dxx.data()[i] += dx.data()[i2]*dx.data()[i2];
			dyy.data()[i] += dy.data()[i2]*dy.data()[i2];
			dxy.data()[i] += dx.data()[i2]*dy.data()[i2];
		}
	// Smooth second moment matrix
	NFilter::recursiveSmoothX(dxx,aRho);
	NFilter::recursiveSmoothY(dxx,aRho);
	NFilter::recursiveSmoothX(dyy,aRho);
	NFilter::recursiveSmoothY(dyy,aRho);
	NFilter::recursiveSmoothX(dxy,aRho);
	NFilter::recursiveSmoothY(dxy,aRho);
	// Compute smallest eigenvalue
	for (int i = 0; i < aSize; i++) {
		float a = dxx.data()[i];
		float b = dxy.data()[i];
		float c = dyy.data()[i];
		float temp = 0.5*(a+c);
		float temp2 = temp*temp+b*b-a*c;
		if (temp2 < 0.0f) aCorners.data()[i] = 0.0f;
		else aCorners.data()[i] = temp-sqrt(temp2);
	}
}

// Code fragment from Pedro Felzenszwalb	---------------
// http://people.cs.uchicago.edu/~pff/dt/
void dt(CVector<float>& f, CVector<float>& d, int n) {
	d.setSize(n);
	int *v = new int[n];
	float *z = new float[n+1];
	int k = 0;
	v[0] = 0;
	z[0] = -10e20;
	z[1] = 10e20;
	for (int q = 1; q <= n-1; q++) {
		float s	= ((f[q]+q*q)-(f(v[k])+v[k]*v[k]))/(2*(q-v[k]));
		while (s <= z[k]) {
			k--;
			s	= ((f[q]+q*q)-(f(v[k])+v[k]*v[k]))/(2*(q-v[k]));
		}
		k++;
		v[k] = q;
		z[k] = s;
		z[k+1] = 10e20;
	}
	k = 0;
	for (int q = 0; q <= n-1; q++) {
		while (z[k+1] < q)
			k++;
		int help = q-v[k];
		d(q) = help*help + f(v[k]);
	}
	delete[] v;
	delete[] z;
}
// ------------------------------------------------------

// euclideanDistanceTransform
void euclideanDistanceTransform(CMatrix<float>& aMatrix) {
	int aXSize = aMatrix.xSize();
	int aYSize = aMatrix.ySize();
	CVector<float> f(NMath::max(aXSize,aYSize));
	// Transform along columns
	for (int x = 0; x < aXSize; x++) {
		for (int y = 0; y < aYSize; y++)
			f(y) = aMatrix(x,y);
		CVector<float> d;
		dt(f,d,aYSize);
		for (int y = 0; y < aYSize; y++)
			aMatrix(x,y) = d(y);
	}
	// Transform along rows
	for (int y = 0; y < aYSize; y++) {
		int aOffset = y*aXSize;
		for (int x = 0; x < aXSize; x++)
			f(x) = aMatrix.data()[x+aOffset];
		CVector<float> d;
		dt(f,d,aXSize);
		for (int x = 0; x < aXSize; x++)
			aMatrix.data()[x+aOffset] = d(x);
	}
}

