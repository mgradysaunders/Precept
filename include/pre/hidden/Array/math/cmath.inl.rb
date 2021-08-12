funcs = [
['abs',         [['T','arrx']]],
['arg',         [['T','arrx']]],
['real',        [['T','arrx']]],
['imag',        [['T','arrx']]],
['conj',        [['T','arrx']]],
['norm',        [['T','arrx']]],
#['min',         [['T','arrx'],['T','arry']]],
#['max',         [['T','arrx'],['T','arry']]],
['fabs',        [['T','arrx']]],
['fma',         [['T','arrx'], ['T','arry'], ['T','arrz']]],
['fmin',        [['T','arrx'], ['T','arry']]],
['fmax',        [['T','arrx'], ['T','arry']]],
['fdim',        [['T','arrx'], ['T','arry']]],
['fmod',        [['T','arrx'], ['T','arry']]],
['remquo',      [['T','arrx'], ['T','arry'], ['int*','arrq']]],
['remainder',   [['T','arrx'], ['T','arry']]],
['nearbyint',   [['T','arrx']]],
['floor',       [['T','arrx']]],
['ceil',        [['T','arrx']]],
['trunc',       [['T','arrx']]],
['round',       [['T','arrx']]],
['rint',        [['T','arrx']]],
['lrint',       [['T','arrx']]],
['llrint',      [['T','arrx']]],
['lround',      [['T','arrx']]],
['llround',     [['T','arrx']]],
['frexp',       [['T','arrx'], ['int*','arrp']]],
['ldexp',       [['T','arrx'], ['int','arrp']]],
['logb',        [['T','arrx']]],
['ilogb',       [['T','arrx']]],
['scalbn',      [['T','arrx'], ['int','arrp']]],
['scalbln',     [['T','arrx'], ['long','arrp']]],
['modf',        [['T','arrx'], ['T*','arrp']]],
['nextafter',   [['T','arrx'], ['T','arry']]],
['nexttoward',  [['T','arrx'], ['long double','arry']]],
['copysign',    [['T','arrx'], ['T','arry']]],
['signbit',     [['T','arrx']]],
['isnan',       [['T','arrx']]],
['isinf',       [['T','arrx']]],
['isfinite',    [['T','arrx']]],
['isnormal',    [['T','arrx']]],
['exp',         [['T','arrx']]],
['log',         [['T','arrx']]],
['exp2',        [['T','arrx']]],
['log2',        [['T','arrx']]],
['log10',       [['T','arrx']]],
['expm1',       [['T','arrx']]],
['log1p',       [['T','arrx']]],
['pow',         [['T','arrx'], ['T','arry']]],
['sqrt',        [['T','arrx']]],
['cbrt',        [['T','arrx']]],
['hypot',       [['T','arrx'], ['T','arry']]],
['erf',         [['T','arrx']]],
['erfc',        [['T','arrx']]],
['lgamma',      [['T','arrx']]],
['tgamma',      [['T','arrx']]],
['sin',         [['T','arrx']]],
['cos',         [['T','arrx']]],
['tan',         [['T','arrx']]],
['asin',        [['T','arrx']]],
['acos',        [['T','arrx']]],
['atan',        [['T','arrx']]],
['atan2',       [['T','arry'], ['T','arrx']]],
['sinh',        [['T','arrx']]],
['cosh',        [['T','arrx']]],
['tanh',        [['T','arrx']]],
['asinh',       [['T','arrx']]],
['acosh',       [['T','arrx']]],
['atanh',       [['T','arrx']]],
['csc',         [['T','arrx']]],
['sec',         [['T','arrx']]],
['cot',         [['T','arrx']]],
['csch',        [['T','arrx']]],
['sech',        [['T','arrx']]],
['coth',        [['T','arrx']]],
['acsc',        [['T','arrx']]],
['asec',        [['T','arrx']]],
['acot',        [['T','arrx']]],
['acsch',       [['T','arrx']]],
['asech',       [['T','arrx']]],
['acoth',       [['T','arrx']]]
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
