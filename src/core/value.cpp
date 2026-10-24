#include "mdbridge/core/value.hpp"

#include <cstdio>
#include <cstdlib>

namespace mdbridge::core {

std::string DateTime::isoformat() const {
    char date[32];
    std::snprintf(date, sizeof date, "%04d-%02d-%02d", year, month, day);
    char time[64];
    if (microsecond != 0) {
        std::snprintf(time, sizeof time, "%02d:%02d:%02d.%06d", hour, minute, second, microsecond);
    } else {
        std::snprintf(time, sizeof time, "%02d:%02d:%02d", hour, minute, second);
    }
    switch (kind) {
    case Kind::LocalDate:
        return date;
    case Kind::LocalTime:
        return time;
    case Kind::LocalDateTime:
        return std::string(date) + "T" + time;
    case Kind::OffsetDateTime:
        break;
    }
    char offset[16];
    const int minutes = std::abs(offsetMinutes);
    std::snprintf(offset, sizeof offset, "%c%02d:%02d", offsetMinutes < 0 ? '-' : '+', minutes / 60,
                  minutes % 60);
    return std::string(date) + "T" + time + offset;
}

const Value* Value::find(std::string_view key) const {
    const Table* table = get<Table>();
    if (!table) {
        return nullptr;
    }
    for (const auto& [k, v] : *table) {
        if (k == key) {
            return &v;
        }
    }
    return nullptr;
}

Value* Value::find(std::string_view key) {
    return const_cast<Value*>(std::as_const(*this).find(key));
}

const Value* Value::path(std::initializer_list<std::string_view> keys) const {
    const Value* node = this;
    for (std::string_view key : keys) {
        node = node->find(key);
        if (!node) {
            return nullptr;
        }
    }
    return node;
}

std::size_t Value::size() const {
    if (const Array* array = get<Array>()) {
        return array->size();
    }
    if (const Table* table = get<Table>()) {
        return table->size();
    }
    return 0;
}

std::string Value::stringOr(std::string_view fallback) const {
    if (const std::string* s = get<std::string>()) {
        return *s;
    }
    return std::string(fallback);
}

std::string_view typeName(Value::Type type) {
    switch (type) {
    case Value::Type::Null:
        return "none";
    case Value::Type::Bool:
        return "boolean";
    case Value::Type::Int:
        return "integer";
    case Value::Type::Float:
        return "float";
    case Value::Type::String:
        return "string";
    case Value::Type::DateTime:
        return "date-time";
    case Value::Type::Array:
        return "array";
    case Value::Type::Table:
        return "table";
    }
    return "value";
}

}  // namespace mdbridge::core
