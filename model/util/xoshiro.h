/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2026 The Kids Research Institute Australia
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

/*
 * RNG is from: http://prng.di.unimi.it/
 */

#ifndef OM_util_xoshiro
#define OM_util_xoshiro

#include <cstdint>
#include <limits>

/// Implementation of Xoshiro256+
class Xoshiro256P {
public:
    typedef uint64_t result_type;
    
    Xoshiro256P(uint64_t a, uint64_t b = 0x220b5e8d, uint64_t c = 0x712a58a2, uint64_t d = 0x712a58a2) {
        this->seed(a, b, c, d);
    }

//     template<class Sseq>
//     explicit Xoshiro256P(Sseq& seq);
    
    template<typename R>
    explicit Xoshiro256P(R& source) {
        seed(source);
    }

    // Disable copying
    Xoshiro256P(const Xoshiro256P&) = delete;
    Xoshiro256P& operator=(const Xoshiro256P&) = delete;

    /// Allow moving, with explicit functions
    Xoshiro256P(Xoshiro256P&& other) {
        memcpy(s, other.s, sizeof s);
    }
    void operator=(Xoshiro256P&& other) {
        memcpy(s, other.s, sizeof s);
    }

//     template<class Sseq> void seed(Sseq& seq);

    void seed(uint64_t a, uint64_t b = 0x220b5e8d, uint64_t c = 0x712a58a2, uint64_t d = 0x712a58a2) {
        s[0] = a;
        s[1] = b;
        s[2] = c;
        s[3] = d;
    }

    template<typename R>
    void seed(R& source) {
        s[0] = source.gen_u64();
        s[1] = source.gen_u64();
        s[2] = source.gen_u64();
        s[3] = source.gen_u64();
    }

    uint64_t operator()();
    uint32_t gen_u32();
    double gen_double();

    friend bool operator==(const Xoshiro256P& lhs, const Xoshiro256P& rhs);
    friend bool operator!=(const Xoshiro256P& lhs, const Xoshiro256P& rhs);

//     template<typename CharT, typename Traits>
//     friend std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os, const Xoshiro256P& rng);
// 
//     template<typename CharT, typename Traits>
//     friend std::basic_istream<CharT, Traits>& operator>>(std::basic_istream<CharT, Traits>& is, Xoshiro256P& rng);

    void binary_checkpoint(ostream& stream) {
        stream.write (reinterpret_cast<char*>(&s[0]), sizeof(s));
    }
    void binary_checkpoint(istream& stream) {
        stream.read (reinterpret_cast<char*>(&s[0]), sizeof(s));
        if (!stream || stream.gcount() != sizeof(s))
            throw runtime_error ("Xoshiro256P::binary_checkpoint: stream read error");
    }

    static constexpr uint64_t min() { return std::numeric_limits<uint64_t>::min(); }
    static constexpr uint64_t max() { return std::numeric_limits<uint64_t>::max(); }

private:
    void generate_block();
    void chacha_core();

    uint64_t s[4];
};

// template<class Sseq> 
// inline Xoshiro256P::Xoshiro256P(Sseq& seq) {
//     seed(seq);
// }
// 
// template<class Sseq>
// inline void Xoshiro256P::seed(Sseq& seq) {
//     seq.generate(s, s + 4);
// }


static inline uint64_t rotl(const uint64_t x, int k) {
    return (x << k) | (x >> (64 - k));
}

inline uint64_t Xoshiro256P::operator()() {
    const uint64_t result = s[0] + s[3];

    const uint64_t t = s[1] << 17;

    s[2] ^= s[0];
    s[3] ^= s[1];
    s[1] ^= s[2];
    s[0] ^= s[3];

    s[2] ^= t;

    s[3] = rotl(s[3], 45);

    return result;
}

inline uint32_t Xoshiro256P::gen_u32() {
    // Since it is recommended not to use the low bits, we take the high part.
    uint64_t x = this->operator()();
    return x >> 32;
}

inline double Xoshiro256P::gen_double() {
    // Doubles support 53-bits of precision (1 bit implied). Use our high bits.
    uint64_t x = this->operator()();
    const double v = 1.1102230246251565e-16; // = 0x1.0p-53
    return (x >> 11) * v;
}


// Implement <random> interface.
inline bool operator==(const Xoshiro256P& lhs, const Xoshiro256P& rhs) {
    for (int i = 0; i < 4; ++i) {
        if (lhs.s[i] != rhs.s[i]) return false;
    }

    return true;
}

inline bool operator!=(const Xoshiro256P& lhs, const Xoshiro256P& rhs) { return !(lhs == rhs); }

// template<typename CharT, typename Traits>
// inline std::basic_ostream<CharT, Traits>& operator<<(std::basic_ostream<CharT, Traits>& os, const Xoshiro256P& rng) {
//     // Serialize.
//     for (int i = 0; i < 4; ++i) os << rng.s[i] << space;
//     os << rng.ctr;
// 
//     // Sestore old state.
//     os.flags(flags);
//     os.fill(fill);
// 
//     return os;
// }
// 
// template<typename CharT, typename Traits>
// inline std::basic_istream<CharT, Traits>& operator>>(std::basic_istream<CharT, Traits>& is, Xoshiro256P& rng) {
//     typedef typename std::basic_istream<CharT, Traits> ::ios_base ios_base;
// 
//     // Save old flags and set ours.
//     auto flags = is.flags();
//     is.flags(ios_base::dec);
// 
//     // Deserialize.
//     for (int i = 0; i < 10; ++i) is >> rng.keysetup[i];
//     is >> rng.ctr;
// 
//     // Restore old flags.
//     is.flags(flags);
// 
//     return is;
// }

#endif
