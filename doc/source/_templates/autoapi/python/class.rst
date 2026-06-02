{% if obj.display %}
   {% set visible_children = obj.children|selectattr("display")|list %}
   {% if is_own_page %}
{{ obj.name }}
{{ "=" * obj.name | length }}

.. py:module:: {{ obj.id }}

{%- set static_attributes = [] -%}
{%- set normal_attributes = [] -%}
{%- for item in visible_children|selectattr("type", "equalto", "attribute")|list -%}
    {%- if item.name is upper -%}
    {%- set static_attributes = static_attributes.append(item) -%}
    {%- else -%}
    {%- set normal_attributes = normal_attributes.append(item) -%}
    {%- endif -%}
{%- endfor -%}

{% if obj.docstring %}


.. autoapi-nested-parse::

   {{ obj.docstring|indent(3) }}
{% endif %}

Constructor
-----------

.. autoapimethod:: {{ obj.id }}.__init__


{% if static_attributes %}
Static attributes
-----------------

.. autoapisummary::
    :nosignatures:

         {% for attribute in static_attributes %}
   {{ attribute.id }}
         {% endfor %}
{% endif %}

{% if normal_attributes %}
Attributes
----------

.. autoapisummary::
    :nosignatures:

         {% for attribute in normal_attributes %}
   {{ attribute.id }}
         {% endfor %}
{% endif %}

{% set visible_methods = visible_children|selectattr("type", "equalto", "method")|list %}
{% if visible_methods %}
Methods
-------

.. autoapisummary::
    :nosignatures:

            {% for method in visible_methods %}
   {{ method.id }}
            {% endfor %}

{% endif %}
{% endif %}
{% endif %}
