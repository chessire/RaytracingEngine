#ifndef __DEEP_RT__
#define __DEEP_RT__

class RaytracingEngine
{
public:
	static RaytracingEngine* Get();

	void Initialize();

private:
	RaytracingEngine();
	~RaytracingEngine();
	void operator=(const RaytracingEngine& deeprt);
};

#endif // __DEEP_RT__
