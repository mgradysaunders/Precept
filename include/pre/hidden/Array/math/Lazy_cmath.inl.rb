funcs = [
['abs',         [['F','arrx']]],
['arg',         [['F','arrx']]],
['real',        [['F','arrx']]],
['imag',        [['F','arrx']]],
['conj',        [['F','arrx']]],
['norm',        [['F','arrx']]],
#['min',         [['F','arrx'],['G','arry']]],
#['max',         [['F','arrx'],['G','arry']]],
['fabs',        [['F','arrx']]],
['fma',         [['F','arrx'], ['G','arry'], ['H','arrz']]],
['fmin',        [['F','arrx'], ['G','arry']]],
['fmax',        [['F','arrx'], ['G','arry']]],
['fdim',        [['F','arrx'], ['G','arry']]],
['fmod',        [['F','arrx'], ['G','arry']]],
#['remquo',      [['F','arrx'], ['G','arry'], ['int*','arrq']]],
['remainder',   [['F','arrx'], ['G','arry']]],
['nearbyint',   [['F','arrx']]],
['floor',       [['F','arrx']]],
['ceil',        [['F','arrx']]],
['trunc',       [['F','arrx']]],
['round',       [['F','arrx']]],
['rint',        [['F','arrx']]],
['lrint',       [['F','arrx']]],
['llrint',      [['F','arrx']]],
['lround',      [['F','arrx']]],
['llround',     [['F','arrx']]],
#['frexp',       [['F','arrx'], ['int*','arrp']]],
#['ldexp',       [['F','arrx'], ['int','arrp']]],
#['logb',        [['F','arrx']]],
#['ilogb',       [['F','arrx']]],
#['scalbn',      [['F','arrx'], ['int','arrp']]],
#['scalbln',     [['F','arrx'], ['long','arrp']]],
#['modf',        [['F','arrx'], ['T*','arrp']]],
['nextafter',   [['F','arrx'], ['G','arry']]],
['nexttoward',  [['F','arrx'], ['G','arry']]],
['copysign',    [['F','arrx'], ['G','arry']]],
#['signbit',     [['F','arrx']]],
#['isnan',       [['F','arrx']]],
#['isinf',       [['F','arrx']]],
#['isfinite',    [['F','arrx']]],
#['isnormal',    [['F','arrx']]],
['exp',         [['F','arrx']]],
['log',         [['F','arrx']]],
['exp2',        [['F','arrx']]],
['log2',        [['F','arrx']]],
['log10',       [['F','arrx']]],
['expm1',       [['F','arrx']]],
['log1p',       [['F','arrx']]],
['pow',         [['F','arrx'], ['G','arry']]],
['sqrt',        [['F','arrx']]],
['cbrt',        [['F','arrx']]],
['hypot',       [['F','arrx'], ['G','arry']]],
['erf',         [['F','arrx']]],
['erfc',        [['F','arrx']]],
['lgamma',      [['F','arrx']]],
['tgamma',      [['F','arrx']]],
['sin',         [['F','arrx']]],
['cos',         [['F','arrx']]],
['tan',         [['F','arrx']]],
['asin',        [['F','arrx']]],
['acos',        [['F','arrx']]],
['atan',        [['F','arrx']]],
['atan2',       [['G','arry'], ['F','arrx']]],
['sinh',        [['F','arrx']]],
['cosh',        [['F','arrx']]],
['tanh',        [['F','arrx']]],
['asinh',       [['F','arrx']]],
['acosh',       [['F','arrx']]],
['atanh',       [['F','arrx']]],
['csc',         [['F','arrx']]],
['sec',         [['F','arrx']]],
['cot',         [['F','arrx']]],
['csch',        [['F','arrx']]],
['sech',        [['F','arrx']]],
['coth',        [['F','arrx']]],
['acsc',        [['F','arrx']]],
['asec',        [['F','arrx']]],
['acot',        [['F','arrx']]],
['acsch',       [['F','arrx']]],
['asech',       [['F','arrx']]],
['acoth',       [['F','arrx']]]
]

puts <<STR
namespace pre {

STR

for func in funcs
    funcname = func[0]
    args1 = []
    args2 = []
    tparams = []
    for arg in func[1]
        args1 << "const LazyArray<#{arg[0]}, R>& #{arg[1]}"
        args2 << "#{arg[1]}(k)"
        tparams << "typename #{arg[0]}"
    end
    args1 = args1.join ", "
    args2 = args2.join ", "
    tparams = tparams.join ", "
    puts <<STR
template <#{tparams}, size_t R>
[[gnu::always_inline]]
inline auto #{funcname}(#{args1}) noexcept
{
    return LazyArray([&](auto k) { return pre::#{funcname}(#{args2}); }, #{func[1][0][1]}.sizes);
}

STR
end

puts <<STR

} // namespace pre

STR
