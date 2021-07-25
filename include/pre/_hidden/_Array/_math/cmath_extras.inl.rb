funcs = [
['min',        [['T','arrx'], ['T','arry']]],
['max',        [['T','arrx'], ['T','arry']]],
['clamp',      [['T','arrx'], ['T','arra'], ['T','arrb']]],
['nthpow',     [['T','arrx'], ['int','arrn']]],
['rotl',       [['T','arrval'], ['unsigned','arrrot']]],
['rotr',       [['T','arrval'], ['unsigned','arrrot']]],
['float_incr', [['T','arrx']]],
['float_decr', [['T','arrx']]],
['fast_floor', [['T','arrx']]],
['fast_ceil' , [['T','arrx']]],
['fast_trunc', [['T','arrx']]],
['fast_round', [['T','arrx']]],
['fract',      [['T','arrx'], ['int*','arri']]],
['sinpi',      [['T','arrx']]],
['cospi',      [['T','arrx']]],
['erfinv',     [['T','arrx']]]
]

puts <<STR
namespace pre {

STR

for func in funcs
    funcname = func[0]
    args1 = []
    args2 = []
    args3 = []
    incrs = []
    decls = []
    for arg in func[1]
        if arg[0][-1] != '*' 
            args1 << "const Array<#{arg[0]}, N...>& #{arg[1]}"
            args2 << if arg[0] == 'T'
               "#{arg[0]}()"
            else
               "std::declval<#{arg[0]}>()"
            end
            args3 << "*itr#{arg[1]}"
            decls << "auto itr#{arg[1]} = #{arg[1]}.begin();"
        else
            args1 << "Array<#{arg[0].slice(0..-2)}, N...>* #{arg[1]}"
            args2 << "std::declval<#{arg[0]}>()"
            args3 << "&(*itr#{arg[1]})"
            decls << "auto itr#{arg[1]} = #{arg[1]}->begin();"
        end
        incrs << "++itr#{arg[1]}"
    end
    args1 = args1.join ", "
    args2 = args2.join ", "
    args3 = args3.join ", "
    incrs << "++itrres"
    incrs = incrs.join ", "
    decls = decls.join "\n"
    restype = "Array<decltype(pre::#{funcname}(#{args2})), N...>"
    puts <<STR
template <typename T, size_t... N>
inline auto #{funcname}(#{args1}) noexcept
{
    #{restype} res;
    #{decls}
    auto itrres = res.begin();
    for (; itrres < res.end(); #{incrs})
        *itrres = pre::#{funcname}(#{args3});
    return res;
}

STR
end

puts <<STR

} // namespace pre

STR
