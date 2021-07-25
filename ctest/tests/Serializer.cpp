#include "../doctest.h"
#include <sstream>
#include <pre/BoundBox>
#include <pre/IdString>
#include <pre/Serializer>

static int Serial_ctor_calls = 0;
static int Serial_dtor_calls = 0;
static int Serial_serialize_calls = 0;

class Serial : virtual public pre::Serializable {
  public:
    Serial() {
        ++Serial_ctor_calls;
    }
    ~Serial() {
        ++Serial_dtor_calls;
    }
};

class SerialBranch final : public Serial {
  public:
    SERIALIZABLE(SerialBranch);
    void serialize(pre::Serializer& serializer) {
        serializer <=> serials;
        ++Serial_serialize_calls;
    }

  public:
    std::map<std::string, pre::RefPtr<Serial>> serials;
};

class SerialLeaf final : public Serial {
  public:
    SERIALIZABLE(SerialLeaf);
    void serialize(pre::Serializer& serializer) {
        serializer <=> str;
        serializer <=> box;
        ++Serial_serialize_calls;
    }

    SerialLeaf() = default;
    SerialLeaf(const pre::IdString& str_, const pre::BoundBox3<double>& box_)
        : str(str_), box(box_) {
    }

  public:
    pre::IdString str;
    pre::BoundBox3<double> box;
};

TEST_CASE("Serializer") {
    std::stringstream ss;
    {
        pre::RefPtr branch1(new SerialBranch());
        pre::RefPtr branch2(new SerialBranch());
        pre::RefPtr leaf1(new SerialLeaf(
                "Leaf1", pre::BoundBox{
                                 pre::Array{-1.0, -2.0, -3.0},
                                 pre::Array{+1.0, +2.0, +3.0}}));
        pre::RefPtr leaf2(new SerialLeaf(
                "Leaf2", pre::BoundBox{
                                 pre::Array{-3.0, -2.0, -1.0},
                                 pre::Array{+4.0, +5.0, +6.0}}));
        pre::RefPtr leaf3(new SerialLeaf(
                "Leaf3", pre::BoundBox{
                                 pre::Array{-7.0, -7.0, -5.0},
                                 pre::Array{+1.0, +1.0, +1.0}}));
        branch1->serials["1"] = leaf1;
        branch1->serials["2"] = leaf2;
        branch2->serials["3"] = leaf3;
        branch2->serials["1"] = leaf1;
        branch2->serials["2"] = branch1;
        pre::StandardSerializer serializer(static_cast<std::ostream&>(ss));
        serializer <=> branch2;
        CHECK(Serial_ctor_calls == Serial_serialize_calls);
    }

    {
        Serial_ctor_calls = 0;
        Serial_dtor_calls = 0;
        pre::RefPtr<Serial> foo(nullptr);
        pre::StandardSerializer serializer(static_cast<std::istream&>(ss));
        serializer <=> foo;
        CHECK(Serial_ctor_calls == 5);
    }

    CHECK(Serial_ctor_calls == Serial_dtor_calls);
}
