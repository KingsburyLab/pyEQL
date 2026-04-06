#define F_MAIN FC_FUNC(f_main, F_MAIN)

extern "C" int F_MAIN();

int main(void)
{
  return F_MAIN();
}
