#ifndef _VECTOR_H_
#define _VECTOR_H_

#include <cstdint>
#include <stdexcept>
#include <utility>
#include <iterator>

/*
    Utility gives std::rel_ops which will fill in relational
        iterator operations so long as operators are provided.
    In any case, ensure that all operations listed in the below website
        are legal for use with iterators.

    http://www.cplusplus.com/reference/iterator/RandomAccessIterator/
 */

using namespace std::rel_ops;

namespace epl
{
    class Invalid_Iterator
    {
        public:
            enum SeverityLevel
            {
                SEVERE,
                MODERATE,
                MILD,
                WARNING
            };

            SeverityLevel level;

            Invalid_Iterator(SeverityLevel level = SEVERE) { this->level = level; }

            virtual const char* what() const
            {
                switch (level)
                {
                    case WARNING:
                        return "Warning";

                    case MILD:
                        return "Mild";

                    case MODERATE:
                        return "Moderate";

                    case SEVERE:
                        return "Severe";

                    default:
                        return "ERROR";
                }
            }
    };

    template <typename T>
    class Vector 
    {
        private:
            // Constant
            static constexpr uint64_t MIN_CAPACITY = 8;

            // Data member
            T* cap; // Start address of the allocated storage
            T* len; // Start address of constructed elements in the vector 
            uint64_t capacity; // Storage size
            uint64_t length; // Vector size (number of elements)
            uint64_t modification; // Version number for vector modification
            uint64_t reallocation; // Version number for memory reallocation

            // Nested class: iterator
            template <bool is_const_iterator = true>
            class Iterator_Base : public std::iterator<std::random_access_iterator_tag, T>
            {
                private:
                    // Type alias
                    using data_structure_pointer_type = typename std::conditional<is_const_iterator, Vector const*, Vector*>::type;
                    using value_reference_type = typename std::conditional<is_const_iterator, T const&, T&>::type;
                    using value_pointer_type = typename std::conditional<is_const_iterator, T const*, T*>::type;

                    // Data member
                    data_structure_pointer_type ds;
                    int64_t offset;
                    uint64_t modification_snapshot;
                    uint64_t reallocation_snapshot;

                public:
                    // Type alias
                    using typename std::iterator<std::random_access_iterator_tag, T>::iterator_category;
                    using typename std::iterator<std::random_access_iterator_tag, T>::value_type;
                    using typename std::iterator<std::random_access_iterator_tag, T>::difference_type;
                    using typename std::iterator<std::random_access_iterator_tag, T>::pointer;
                    using typename std::iterator<std::random_access_iterator_tag, T>::reference;

                    // Copy-constructible, copy-assignable and destructible
                    Iterator_Base() : ds{ nullptr },
                                      offset{ 0 },
                                      modification_snapshot{ 0 },
                                      reallocation_snapshot{ 0 } {}

                    Iterator_Base(data_structure_pointer_type ds, int64_t offset, uint64_t modification, uint64_t reallocation) : ds{ ds },
                                                                                                                                  offset{ offset },
                                                                                                                                  modification_snapshot{ modification },
                                                                                                                                  reallocation_snapshot{ reallocation } {}

                    Iterator_Base(Iterator_Base const& that) { Iterator_Base::copy(that); }

                    Iterator_Base& operator=(Iterator_Base const& that)
                    {
                        if (this != &that)
                        {
                            copy(that);
                        }

                        return *this;
                    }

                    ~Iterator_Base() {}

                    /*
                        An iterator can be converted (without warning or
                            type casts) to const_iterator.
                     */                                   
                    operator Iterator_Base<true>(void)
                    {
                        return Iterator_Base<true>(ds, offset, modification_snapshot, reallocation_snapshot);
                    }

                    // Can be incremented
                    Iterator_Base& operator++(void) // Pre-increment operator, e.g. ++it 
                    {
                        assert_valid();
                        ++offset;
                        return *this;
                    }                

                    Iterator_Base operator++(int) // Post-increment operator, e.g. it++
                    {
                        Iterator_Base t{ *this };
                        operator++();
                        return t;
                    }            

                    // Supports equality/inequality comparisons
                    bool operator==(Iterator_Base const& that) const
                    {
                        assert_valid();
                        that.assert_valid();
                        return ds == that.ds && offset == that.offset;
                    }                                         

                    // Can be dereferenced
                    value_reference_type operator*(void) const
                    {
                        assert_valid();
                        return (*ds)[offset];
                    }

                    value_pointer_type operator->(void) const
                    {
                        return &(operator*());
                    }

                    // Can be decremented
                    Iterator_Base& operator--(void) // Pre-decrement operator, e.g. --it
                    {
                        assert_valid();
                        --offset;
                        return *this;
                    }

                    Iterator_Base operator--(int) // Post-decrement operator, e.g. it++
                    {
                        Iterator_Base t{ *this };
                        operator--();
                        return t;
                    }

                    // Supports compound assignment operations += and -=
                    Iterator_Base& operator+=(int64_t k)
                    {
                        assert_valid();
                        offset += k;
                        return *this;
                    }

                    Iterator_Base& operator-=(int64_t k)
                    {
                        operator+=(-k);
                        return *this;
                    }

                    // Supports arithmetic operators + and -
                    Iterator_Base operator+(int64_t k) const
                    {
                        Iterator_Base t{ *this };
                        return t += k;
                    }

                    Iterator_Base operator-(int64_t k) const
                    {
                        Iterator_Base t{ *this };
                        return t -= k;
                    }

                    difference_type operator-(Iterator_Base const& that) const
                    {
                        assert_valid();
                        that.assert_valid();
                        return offset - that.offset;
                    }

                    // Supports inequality comparisons (<, >, <= and >=) between iterators 
                    bool operator<(Iterator_Base const& that) const
                    {
                        return *this - that < 0;
                    }

                    // Supports offset dereference operator ([])
                    value_reference_type operator[](uint64_t k) const
                    {
                        return *(*this + k);
                    }

                private:
                    void copy(Iterator_Base const& that)
                    {
                        ds = that.ds;
                        offset = that.offset;
                        modification_snapshot = that.modification_snapshot;
                        reallocation_snapshot = that.reallocation_snapshot;
                    }

                    void assert_valid(void) const
                    {
                        bool is_modified = (modification_snapshot != ds->modification);
                        bool is_reallocated = (reallocation_snapshot != ds->reallocation);
                        bool is_inbound = (offset >= 0 && offset < ds->length);

                        if (is_modified || is_reallocated)
                        {
                            if (!is_inbound)
                                throw Invalid_Iterator(Invalid_Iterator::SEVERE);
                            else if (is_reallocated)
                                throw Invalid_Iterator(Invalid_Iterator::MODERATE);
                            else 
                                throw Invalid_Iterator(Invalid_Iterator::MILD);
                        }
                    }
            };

        public:
            using const_iterator = Iterator_Base<true>;
            using iterator = Iterator_Base<false>;

            Vector(void)
            {
                init(MIN_CAPACITY);
            }

            explicit Vector(uint64_t n)
            {
                if (n)
                {
                    init(n);
                    length = n;

                    for (uint64_t i = 0; i < n; i += 1)
                    {
                        new (len + i) T{};
                    }
                }
                else 
                {
                    Vector();
                }
            }

            Vector(Vector const& that)
            {
                copy(that);
            }

            Vector(Vector&& that)
            {
                move(std::move(that));
            }

            template <typename IT>
            Vector(IT b, IT e) : Vector(b, e, typename std::iterator_traits<IT>::iterator_category{}) {}

            template <typename IT>
            Vector(IT b, IT e, std::random_access_iterator_tag)
            {
                init(e - b);
                length = 0;

                while (b != e)
                {
                    new (len + length) T{ *b };
                    ++b;
                    ++length;
                }
            }

            template <typename IT>
            Vector(IT b, IT e, std::input_iterator_tag) : Vector() 
            {
                while (b != e)
                {
                    push_back(*b);
                    ++b;
                }
            }

            Vector(std::initializer_list<T> list) : Vector(list.begin(), list.end()) {}

            Vector& operator=(Vector const& that)
            {
                if (this != &that)
                {
                    destroy();
                    copy(that);
                }

                return *this;
            }

            Vector& operator=(Vector&& that)
            {
                if (this != &that)
                {
                    destroy();
                    move(std::move(that));
                }

                return *this;
            }

            ~Vector() 
            {
                destroy();
            }

            // Member function group: Capacity 
            uint64_t size(void) const
            {
                return length;
            }

            // Member function group: Element Access 
            T const& operator[](uint64_t n) const
            {
                if (n >= length)
                {
                    throw std::out_of_range{ "index is out of range" };
                }

                return len[n];
            }

            T& operator[](uint64_t n)
            {
                return const_cast<T&>(static_cast<Vector const&>(*this)[n]);
            }

            // Member function group: Modifiers
            void pop_back(void)
            {
                if (!length)
                {
                    throw std::out_of_range{ "fatal error: popping from an empty vector" };
                }

                (len + length - 1) -> ~T();
                length -= 1;
                modification += 1;
            }

            void pop_front(void)
            {
                if (!length)
                {
                    throw std::out_of_range{ "fatal error: popping from an empty vector" };
                }

                len -> ~T();
                len++;
                length -= 1;
                modification += 1;
            }

            /*
                Since ensure_back_capacity will destroy original vector,
                    if it is the element of original vector that is to be push_back
                    a segmentation fault will occur.
                In-turn, a temp variable is constructed.

                Example:
                    vector<vector<Foo>> x(8);
                    x[0].push_front(Foo());
                    x.push_back(x[0]);
             */
            void push_back(T const& elem)
            {
                T temp{ elem };
                ensure_back_capacity(1);
                new (len + length) T{ std::move(temp) };
                length += 1;
                modification += 1;
            }

            void push_back(T&& elem)
            {
                T temp{ std::move(elem) };
                ensure_back_capacity(1);
                new (len + length) T{ std::move(temp) };
                length += 1;
                modification += 1;
            }

            void push_front(T const& elem)
            {
                T temp{ elem };
                ensure_front_capacity(1);
                new (len - 1) T{ std::move(temp) };
                len--;
                length += 1;
                modification += 1;
            }

            void push_front(T&& elem)
            {
                T temp{ std::move(elem) };
                ensure_front_capacity(1);
                new (len - 1) T{ std::move(temp) };
                len--;
                length += 1;
                modification += 1;
            }

            template <typename... Args>
            void emplace_back(Args... args)
            {
                ensure_back_capacity(1);
                new (len + length) T{ std::forward<Args>(args)... };
                length += 1;
                modification += 1;
            }

            // Member function group: Iterators 
            const_iterator begin(void) const
            {
                return const_iterator(this, 0, modification, reallocation);
            }

            const_iterator end(void) const
            {
                return const_iterator(this, length, modification, reallocation);
            }

            iterator begin(void)
            {
                return iterator(this, 0, modification, reallocation);
            }

            iterator end(void)
            {
                return iterator(this, length, modification, reallocation);
            }

        private:
            void init(uint64_t n)
            {
                cap = static_cast<T*>(::operator new(n * sizeof(T)));
                len = cap;
                capacity = n;
                length = 0;
                modification = 0;
                reallocation = 0;
            }

            void copy(Vector const& that)
            {
                capacity = that.capacity;
                length = that.length;
                cap = static_cast<T*>(::operator new(capacity * sizeof(T)));
                len = cap + (that.len - that.cap);

                for (uint64_t i = 0; i < length; i += 1)
                {
                    new (len + i) T{ that[i] };
                }

                reallocation += 1;
            }

            void move(Vector&& that)
            {
                capacity = that.capacity;
                length = that.length;
                cap = that.cap;
                len = that.len;
                that.cap = nullptr;
                that.len = nullptr;

                /*
                    The storage cap points to is different,
                        it can be considered as a kind of "reallocation".
                 */
                reallocation += 1;
                that.reallocation += 1;
            }

            void destroy(void)
            {
                if (len)
                {
                    for (uint64_t i = 0; i < length; i += 1)
                    {
                        (len + i) -> ~T();
                    }
                }

                if (cap)
                {
                    ::operator delete(cap);
                }
            }

            void ensure_back_capacity(uint64_t capacity_required)
            {
                uint64_t back_capacity = (cap + capacity) - (len + length);

                if (back_capacity >= capacity_required)
                    return;

                uint64_t new_capacity = capacity;

                while (back_capacity < capacity_required)
                {
                    new_capacity += capacity;
                    back_capacity += capacity;
                }

                T* new_cap = static_cast<T*>(::operator new(new_capacity * sizeof(T)));
                T* new_len = new_cap + (len - cap);

                for (uint64_t i = 0; i < length; i += 1)
                {
                    new (new_len + i) T{ std::move(len[i]) };
                    (len + i) -> ~T();
                }

                ::operator delete(cap);

                cap = new_cap;
                len = new_len;
                capacity = new_capacity;
                reallocation += 1;
            }

            void ensure_front_capacity(uint64_t capacity_required)
            {
                uint64_t front_capacity = len - cap;

                if (front_capacity >= capacity_required)
                    return;

                uint64_t new_capacity = capacity;

                while (front_capacity < capacity_required)
                {
                    new_capacity += capacity;
                    front_capacity += capacity;
                }

                T* new_cap = static_cast<T*>(::operator new(new_capacity * sizeof(T)));
                T* new_len = new_cap + new_capacity - length - ((cap + capacity) - (len + length));

                for (uint64_t i = 0; i < length; i += 1)
                {
                    new (new_len + i) T{ std::move(len[i]) };
                    (len + i) -> ~T();
                }

                ::operator delete(cap);

                cap = new_cap;
                len = new_len;
                capacity = new_capacity;
                reallocation += 1;
            }
    };
} // End namespace epl

#endif