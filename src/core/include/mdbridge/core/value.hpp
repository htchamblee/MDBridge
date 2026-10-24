#pragma once
// Plain data of MDATK's TOML files: the tree of tables, arrays and scalars that MDATK's
// `tomlio` reads and writes (Python dicts, lists, str, int, float, bool, datetime).

#include <cstdint>
#include <initializer_list>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace mdbridge::core {

// A TOML date, time or date-time. Fractional seconds are kept to the microsecond, as Python's
// tomllib does (it truncates finer digits).
struct DateTime {
    enum class Kind { LocalDate, LocalTime, LocalDateTime, OffsetDateTime };

    Kind kind = Kind::LocalDateTime;
    int year = 0, month = 0, day = 0;
    int hour = 0, minute = 0, second = 0, microsecond = 0;
    int offsetMinutes = 0;  // OffsetDateTime only

    // The text Python's isoformat() gives, e.g. "2026-10-15T14:32:05" or "...05.250000+02:00".
    std::string isoformat() const;

    bool operator==(const DateTime&) const = default;
};

class Value;
using Array = std::vector<Value>;
// Keys keep the order of the file (Python dicts are ordered).
using Table = std::vector<std::pair<std::string, Value>>;

class Value {
public:
    enum class Type { Null, Bool, Int, Float, String, DateTime, Array, Table };

    Value() = default;
    Value(bool v) : m_data(v) {}
    Value(int v) : m_data(std::int64_t{v}) {}
    Value(std::int64_t v) : m_data(v) {}
    Value(double v) : m_data(v) {}
    Value(const char* v) : m_data(std::string(v)) {}
    Value(std::string v) : m_data(std::move(v)) {}
    Value(DateTime v) : m_data(v) {}
    Value(Array v) : m_data(std::move(v)) {}
    Value(Table v) : m_data(std::move(v)) {}

    Type type() const { return static_cast<Type>(m_data.index()); }
    bool isNull() const { return type() == Type::Null; }
    bool isTable() const { return type() == Type::Table; }
    bool isArray() const { return type() == Type::Array; }
    bool isString() const { return type() == Type::String; }

    // The stored value, or nullptr if this value has another type.
    template <class T>
    const T* get() const { return std::get_if<T>(&m_data); }
    template <class T>
    T* get() { return std::get_if<T>(&m_data); }

    // Table entry `key`, or nullptr (also when this is not a table).
    const Value* find(std::string_view key) const;
    Value* find(std::string_view key);
    // Nested table entry, e.g. path({"inputs", "topology"}).
    const Value* path(std::initializer_list<std::string_view> keys) const;

    // Number of array elements or table entries (0 for scalars).
    std::size_t size() const;

    // A string scalar, or `fallback` for anything else.
    std::string stringOr(std::string_view fallback = {}) const;

    bool operator==(const Value&) const = default;

private:
    std::variant<std::monostate, bool, std::int64_t, double, std::string, DateTime, Array, Table> m_data;
};

// Name of a type as MDATK's messages use it ("table", "array", "string", ...).
std::string_view typeName(Value::Type type);

}  // namespace mdbridge::core
