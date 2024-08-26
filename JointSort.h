#include <cstddef>
#include <iterator>
#include <algorithm>

#ifdef _WIN32
using difference_type = std::ptrdiff_t;
#else
#include <sys/types.h>
using difference_type = ssize_t;
#endif

template <class Key, class Val>
struct JSIterator {
    Key *key;
    Val *val;

    // This is the value type being iterated over.
    struct Pair {
        Key key;
        Val val;

        // If you provide a comparator that takes key, we can convert
        operator Key () const { return key; }

        // If Key has an operator defined, this is used.
        bool operator < (const Pair &other) const { return key < other.key; }
    };

    using value_type = Pair;

    // This represents a reference to Pair. The reason we create this as a
    // new class type, rather than simply using Pair& is that we get to
    struct reference {
        Key *key;
        Val *val;

        reference& operator =(const value_type& v) {
            *key = v.key;
            *val = v.val;
            return *this;
        }

        reference& operator =(const reference& other) {
            *key = *other.key;
            *val = *other.val;
            return *this;
        }

        operator value_type() const {
            return { *key, *val };
        }

        operator Key () const { return *key; }

        bool operator < (const reference& other) const { return *key < *other.key; }
        bool operator < (const value_type& v) const { return *key < v.key; }

        // Makes swap available for ADL lookup.
        friend void swap(reference a, reference b) {
            std::swap(*a.key, *b.key);
            std::swap(*a.val, *b.val);
        };
    };

    using difference_type = difference_type;
    using pointer = reference;
    using iterator_category = std::random_access_iterator_tag;

    inline JSIterator& operator ++() {
        ++key;
        ++val;
        return *this;
    }

    inline JSIterator& operator --() {
        --key;
        --val;
        return *this;
    }

    inline JSIterator& operator +=(difference_type n) {
        key += n;
        val += n;
        return *this;
    }

    inline JSIterator& operator -=(difference_type n) {
        key -= n;
        val -= n;
        return *this;
    }

    inline bool operator != (const JSIterator& other) const {
        return key != other.key;
    }

    inline bool operator == (const JSIterator& other) const {
        return key == other.key;
    }

    inline bool operator < (const JSIterator& other) const {
        return key < other.key;
    }

    inline bool operator > (const JSIterator& other) const {
        return key > other.key;
    }

    inline bool operator <= (const JSIterator& other) const {
        return key <= other.key;
    }

    inline bool operator >= (const JSIterator& other) const {
        return key >= other.key;
    }

    JSIterator operator +(difference_type n) const {
        return {key + n, val + n};
    }

    JSIterator operator -(difference_type n) const {
        return {key - n, val - n};
    }

    difference_type operator -(const JSIterator& other) const {
        return key - other.key;
    }

    reference operator *() const { return {key, val}; }

private:
    void suppressIncorrectGCCUnusedTypeWarnings_() {
        (void) sizeof(pointer);
        (void) sizeof(iterator_category);
    }
};
