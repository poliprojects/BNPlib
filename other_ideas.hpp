// Argument unpacker:
template <size_t N>
class unpack_caller{
  private:
      template <typename FuncType, size_t... I>
      void call(FuncType &f, std::vector<data_t> &args, indices<I...>){
          f(args[I]...);
      }
  public:
      template <typename FuncType>
      void operator()(FuncType &f, std::vector<data_t> &args){
          call(f, args, BuildIndices<N>{});
      }
};


// Variadic templates:
template<typename... Arguments> class VariadicTemplate; or
template<typename T, typename... Arguments> class VariadicTemplate; or
template<typename... Arguments> void SampleFunction(Arguments..., ...)
