#ifndef __DEEP_RT__
#define __DEEP_RT__

class DeepRT
{
public:
	static DeepRT* Get();

	void Initialize();

private:
	DeepRT();
	~DeepRT();
	void operator=(const DeepRT& deeprt);
};

#endif // __DEEP_RT__
