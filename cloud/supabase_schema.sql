-- Optional Alternative Backend: Supabase (Free Tier PostgreSQL)
-- Paste into Supabase SQL Editor -> Run.
-- Row-Level Security (RLS) allows anonymous INSERT ONLY (participants cannot read other people's rows).

create table if not exists public.trpm8_events (
  event_id text primary key,
  received_at timestamptz not null default now(),
  record_type text not null check (record_type in ('checkin', 'trial')),
  participant_id text not null,
  recorded_at timestamptz not null,
  payload jsonb not null
);

create index if not exists idx_trpm8_events_type_time
  on public.trpm8_events (record_type, recorded_at desc);

create index if not exists idx_trpm8_events_participant
  on public.trpm8_events (participant_id);

alter table public.trpm8_events enable row level security;

-- Allow anonymous PWA clients to INSERT new anonymized observations only
create policy "anon_insert_only"
  on public.trpm8_events
  for insert
  to anon
  with check (
    record_type in ('checkin', 'trial')
    and length(participant_id) between 3 and 40
  );

-- Only authenticated project owner (service_role / SQL editor / authenticated researcher) can SELECT
create policy "owner_read_only"
  on public.trpm8_events
  for select
  to authenticated
  using (true);
